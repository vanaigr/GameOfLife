#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "Misc.h"
#include "Grid.h"
#include <vector>

#include <cassert> 

#include"Timer.h"
#include"AutoTimer.h"
#include"MedianCounter.h"

#include <intrin.h>

#include "data.h"

constexpr size_t sizeOfBatch = 32;// sizeof(__m128i); //dependent on CellsGrid fullSize


bool fieldCell::isAlive(const FieldCell cell) {
	return ((misc::to_underlying(cell)) & 0b1);
}
bool fieldCell::isDead(const FieldCell cell) {
	return !isAlive(cell);
}

FieldCell fieldCell::nextGeneration(const FieldCell cell, const uint32_t aliveNeighboursCount) {
	if (fieldCell::isDead(cell) && aliveNeighboursCount == 3) {
		return FieldCell::ALIVE;
	}
	else if (fieldCell::isAlive(cell) && (aliveNeighboursCount < 2 || aliveNeighboursCount > 3)) {
		return FieldCell::DEAD;
	}
	return cell;
}

static const size_t sizeOfUint32 = sizeof uint32_t;

 class Field::FieldPimpl {
 public:
	 struct current {};
	 struct buffer  {};
public:
	const uint32_t width;
	const uint32_t height;
	const uint32_t size;
private:
	const uint32_t startPadding_uint;
	const uint32_t size_uint;
	const uint32_t endPadding_uint;
	const uint32_t fullSize_uint;
	std::unique_ptr<uint32_t[]> current_; 
	std::unique_ptr<uint32_t[]> buffer_;

public:
	FieldPimpl(const uint32_t gridWidth, const uint32_t gridHeight) :
		width(gridWidth), height(gridHeight), size(gridWidth* gridHeight),
		startPadding_uint{ misc::roundUpIntTo(gridWidth, 32) / 32 }, /*
			top row duplicate before first row 
			+ extra bytes for 
				safety and 
				alignment of current row read in 'void threadUpdateGrid(std::unique_ptr<Field::GridData>&)'
		*/
		size_uint{ misc::roundUpIntTo(gridWidth * gridHeight, 32) / 32 },
		endPadding_uint{ misc::roundUpIntTo(misc::roundUpIntTo(gridWidth, 32) + sizeOfBatch * 8 , 32) / 32 }, /*
			bottom row duplicate after last row 
			+ extra bytes for 
				memcopy/vector instructions and 
				left/right boundary cells update 
			all in 'void threadUpdateGrid(std::unique_ptr<Field::GridData>&)'
		*/
		fullSize_uint{ (misc::roundUpIntTo(gridWidth, 32) / 32) + (misc::roundUpIntTo(gridWidth * gridHeight, 32 * 4) / 32 / 4) + (misc::roundUpIntTo(misc::roundUpIntTo(gridWidth, 32) + sizeOfBatch * 8 , 32) / 32) },
		current_{ new uint32_t[(misc::roundUpIntTo(gridWidth, 32) / 32) + (misc::roundUpIntTo(gridWidth * gridHeight, 32) / 32) + (misc::roundUpIntTo(misc::roundUpIntTo(gridWidth, 32) + sizeOfBatch * 8 , 32) / 32)]{ } },
		buffer_ { new uint32_t[(misc::roundUpIntTo(gridWidth, 32) / 32) + (misc::roundUpIntTo(gridWidth * gridHeight, 32) / 32) + (misc::roundUpIntTo(misc::roundUpIntTo(gridWidth, 32) + sizeOfBatch * 8 , 32) / 32)]{ } }
	{ };
	~FieldPimpl() = default;

	template<typename gridType = current>
	void prepareNext() {
		current_.swap(buffer_);
		//std::memcpy(&current_[0], &(current_[size_uint]), sizeOfUint32 * startPadding_uint); //start pading
		if (startPadding_uint * 32 == width) {
			std::memcpy(&current_[0], &(current_[size_uint]), sizeOfUint32 * startPadding_uint);
		}
		else {
			for (uint32_t i = 0; i < startPadding_uint; i++) {
				current_[startPadding_uint - i] = getCellsu((size_uint + i) * 32);
			}
		}
		if (size % 32 == 0) {
			std::memcpy(&current_[size_uint + startPadding_uint], &(current_[startPadding_uint]), sizeOfUint32 * endPadding_uint);
		}
		else {
			const uint32_t fieldEnd = current_[size_uint + startPadding_uint - 1];
			const uint32_t fieldElements = size % 32;
			const uint32_t paddingElements = 32 - fieldElements;

			const uint32_t end = (fieldEnd & (0xffffffff >> paddingElements));
			current_[size + startPadding_uint - 1] = end;

			for (uint32_t i = 0; i < endPadding_uint; i += 1) {
				const uint32_t cells = current_[startPadding_uint + i];

				current_[size_uint + startPadding_uint + i - 1] += (current_[startPadding_uint] << fieldElements);
				current_[size_uint + startPadding_uint + i] += (current_[startPadding_uint] >> paddingElements);
			}
		}
	}

	template<typename gridType = current>
	void fill(const FieldCell cell) {
		auto* grid{ getGrid<gridType>() };
		const auto val = ~0u * fieldCell::isAlive(cell);
		std::fill(&grid[0], &grid[0] + fullSize_uint, val);
	}

	template<typename gridType = current>
	bool isCellAlive(const int32_t index) const {
		auto* grid{ getGrid<gridType>() };
		const auto index_uint = index / 32;
		const auto shift = misc::mod(index, 32);
		return (grid[index_uint + startPadding_uint] >> shift) & 0b1;
	}

	template<typename gridType = current>
	FieldCell cellAt(const int32_t index) const {
		return isCellAlive<gridType>(index) ? FieldCell::ALIVE : FieldCell::DEAD;
	}

	template<typename gridType = current>
	void setCellAt(const int32_t index, const FieldCell cell) const {
		auto* grid{ getGrid<gridType>() };
		const auto index_uint = index / 32;
		const auto shift = misc::mod(index, 32);
		auto& cur{ grid[index_uint + startPadding_uint] };

		cur = (cur & ~(0b1 << shift)) | (fieldCell::isAlive(cell) << shift);
	}/* tested */

	template<typename gridType = current>
	uint32_t& getCellsInt(int32_t index_int) {
		return getGrid<gridType>()[index_int + startPadding_uint];
	}

	template<typename gridType = current>
	constexpr inline uint32_t getCellsu(int32_t index) {
		auto* grid{ getGrid<gridType>() };
		const int32_t index_int = index / 32;
		if (misc::mod(index, 32) == 0) return getCellsInt<gridType>(index_int);
		const auto pt1 = grid[index_int + startPadding_uint];
		const auto pt2 = grid[index_int + 1 + startPadding_uint];

		const uint32_t offset1 = misc::mod(index, 32);
		const uint32_t offset2 = 32 - offset1;
		return (pt1 >> offset1) | (pt2 << offset2);
	}/* tested */

	template<typename gridType = current>
	constexpr inline void setCellsu(int32_t index, uint32_t cells) {
		auto* grid{ getGrid<gridType>() };
		const int32_t index_int = index / 32;
		if (misc::mod(index, 32) == 0) {
			grid[index_int + startPadding_uint] = cells;
			return;
		}

		const uint32_t offset1 = misc::mod(index, 32);
		const uint32_t offset2 = 32 - offset1;

		grid[index_int + startPadding_uint] = (grid[index_int + startPadding_uint] & ( ~0u >> offset2)) | (cells << offset1);
		grid[index_int + startPadding_uint + 1] = (grid[index_int + startPadding_uint + 1] & (~0u << offset1)) | (cells >> offset2);
	}

 private:
	 template<typename gridType>
	 inline uint32_t* getGrid() const;

	 template<>
	 inline uint32_t* getGrid<current>() const {
		 return this->current_.get();
	 }

	 template<>
	 inline uint32_t* getGrid<buffer>() const {
		 return this->buffer_.get();
	 }
};

 struct Field::GridData {
 private: static const uint32_t samples = 100;
 public:
	 UMedianCounter gridUpdate{ samples }, waiting{ samples }, bufferSend{ samples };
	 unsigned int grid__iteration = 0;

	 uint32_t index;

	 std::unique_ptr<FieldPimpl>& grid;
	 std::atomic_bool& interrupt_flag;
	 std::atomic_bool& gpuBufferLock_flag;
	 uint32_t startBatch;
	 uint32_t endBatch;

	 GLFWwindow* offscreen_context;
	 GLuint bufferP;
	 bool& isWriteOffset;
	 uint32_t offset;

	 uint32_t currentWriteOffset_bytes() {
		 return offset * isWriteOffset;
	 }

	 void generationUpdated() {
		 grid__iteration++;
		 if (grid__iteration % (samples + index) == 0)
			 std::cout << "grid task " << index << ": "
			 << gridUpdate.median() << ' '
			 << waiting.median() << ' '
			 << bufferSend.median() << std::endl;
	 }

 public:
	 GridData(
		 uint32_t index_,
		 std::unique_ptr<FieldPimpl>& grid_,
		 std::atomic_bool& interrupt_flag_,
		 uint32_t startBatch_,
		 uint32_t endBatch_,
		 GLuint bufferP_,
		 bool& isOffset_,
		 uint32_t offset_,
		 std::atomic_bool& gpuBufferLock_flag_,
		 GLFWwindow* offscreen_context_
	 ) :
		 index(index_),
		 grid(grid_),
		 interrupt_flag(interrupt_flag_),
		 gpuBufferLock_flag(gpuBufferLock_flag_),
		 startBatch(startBatch_),
		 endBatch  (endBatch_),
		 offscreen_context(offscreen_context_),
		 bufferP(bufferP_),
		 isWriteOffset(isOffset_),
		 offset(offset_)
	 {}
 };

 static FieldCell updatedCell(const int32_t index, const std::unique_ptr<Field::FieldPimpl>& cellsGrid) {
	 //const bool grid = cellsGrid->cellAt(0);
	 const auto width = cellsGrid->width;
	 const auto height = cellsGrid->height;
	 bool cell = cellsGrid->isCellAlive(index);

	 uint32_t aliveNeighbours = 0;

	 for (int yo = -1; yo <= 1; yo++) {
		 for (int xo = -1; xo <= 1; xo++) {
			 if (xo != 0 || yo != 0) {
				 int x = ((index + xo) + width) % width,
					 y = (((index / width) + yo) + height) % height;
				 int offsetedIndex = x + width * y;
				 if (cellsGrid->isCellAlive(offsetedIndex))
					 aliveNeighbours++;
			 }
		 }
	 }

	 return fieldCell::nextGeneration(cell ? FieldCell::ALIVE : FieldCell::DEAD, aliveNeighbours);
 }

 static void threadUpdateGrid(std::unique_ptr<Field::GridData>& data) {
	 Timer<> t{};
	 auto& grid = data->grid;
	 const int32_t width = static_cast<int32_t>(grid->width);
	 const int32_t height = static_cast<int32_t>(grid->height);
	 const int lastElement = width - 1;

	 GLFWwindow* window = data->offscreen_context;
	 if (window != glfwGetCurrentContext())glfwMakeContextCurrent(window);

	 const auto setBufferCellAt = [&grid](uint32_t index, FieldCell cell) -> void { grid->setCellAt<Field::FieldPimpl::buffer>(index, cell); };

	 struct remainder {
		 uint16_t cellsCols;
		 uint8_t curCell;
	 };

	 const auto calcNewGen_batch = [&grid](
		 const remainder previousRemainder,
		 uint32_t index_batch,
		 remainder& currentRemainder_out
		 ) -> uint32_t {
			 static_assert(sizeOfBatch > 1, "sizeOfBatch = 0, 1 doesnt make sense");
			 const auto width = grid->width;

			 const auto isCellAlive = [&grid](uint32_t index) -> bool { return grid->isCellAlive(index); };
			 const auto put32At4Bits = [](const uint32_t number) -> __m128i {
				 const __m128i indecesLower = _mm_slli_si128(_mm_set1_epi8(1), 8); // lower->higher: 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 
				 const __m128i indecesHigher = _mm_add_epi8(indecesLower, _mm_set1_epi8(2));

				 const __m128i num = _mm_set1_epi32(number);

				 const auto extract4BitsAndMask = [](
					 const __m128i num,
					 const  uint8_t bi1, const  uint8_t bi2, const  uint8_t bi3, const uint8_t bi4, const uint8_t bitMask, const __m128i indeces
					 ) -> __m128i {
					 const auto high_low_bytesToLontTo128 = [](const uint8_t b1, const uint8_t b2, const  uint8_t b3, const uint8_t b4,
						 const uint8_t b5, const uint8_t b6, const  uint8_t b7, const uint8_t b8) -> __m128i {
						 return _mm_set1_epi64x(
							 (uint64_t(b1) << 8 * 0)
							 | (uint64_t(b2) << 8 * 1)
							 | (uint64_t(b3) << 8 * 2)
							 | (uint64_t(b4) << 8 * 3)
							 | (uint64_t(b5) << 8 * 4)
							 | (uint64_t(b6) << 8 * 5)
							 | (uint64_t(b7) << 8 * 6)
							 | (uint64_t(b8) << 8 * 7));
					 };
					 const __m128i mask = high_low_bytesToLontTo128(1, 2, 4, 8, 16, 32, 64, 128);
					 return
						 _mm_and_si128(
							 _mm_cmpeq_epi8(
								 _mm_and_si128(
									 _mm_shuffle_epi8(
										 num,
										 indeces
									 ),
									 mask
								 ),
								 mask
							 ),
							 _mm_set1_epi8(bitMask));
				 };


				 const __m128i resultlower4 = extract4BitsAndMask(num, 0, 2, 4, 6, 0b1, indecesLower);
				 const __m128i resultHigher4 = extract4BitsAndMask(num, 1, 3, 5, 7, 0b1 << 4, indecesHigher);

				 return _mm_or_si128(resultlower4, resultHigher4);
			 };

			 const auto shiftLeft1Data = [](const __m128i data) -> __m128i {
				 const auto carry =
					 _mm_and_si128(
						 _mm_srli_si128(_mm_slli_epi16(data, 4), 15), _mm_set1_epi8(0b11110000u)
					 ); static_assert(
						 "_mm_srli_epi8(_mm_bslli_si128(sum, 15), 4)",
						 "lack of _mm_srli_epi8() !");
				 return _mm_add_epi8(
					 _mm_slli_si128(data, 1),
					 carry
				 );
			 };

			 const auto bits4To32 = [](const __m128i bits4) -> uint32_t {
				 const __m128i maskLower4 = _mm_set1_epi8(0b00000001u);
				 const __m128i maskHigher4 = _mm_set1_epi8(0b00010000u);

				 uint32_t lower16{ static_cast<std::make_unsigned<int>::type>(_mm_movemask_epi8(_mm_cmpeq_epi8(_mm_and_si128(bits4, maskLower4), maskLower4))) };
				 uint32_t higher16{ static_cast<std::make_unsigned<int>::type>(_mm_movemask_epi8(_mm_cmpeq_epi8(_mm_and_si128(bits4, maskHigher4), maskHigher4))) };

				 return lower16 | (higher16 << 16);
			 };

			 static_assert(sizeof FieldCell == 1, "size of field cell must be 1 byte");
			 const uint32_t topCellRow_uint_{ grid->getCellsu(index_batch * sizeOfBatch - width) };
			 const uint32_t curCellRow_uint_{ grid->getCellsInt(index_batch) };
			 const uint32_t botCellRow_uint_{ grid->getCellsu(index_batch * sizeOfBatch + width) };

			 //if ((index % sizeOfUint32) != 0) std::cout << "(sick!)" << index << std::endl;

			 const __m128i topCellRow_ = put32At4Bits(topCellRow_uint_);
			 const __m128i curCellRow_ = put32At4Bits(curCellRow_uint_);
			 const __m128i botCellRow_ = put32At4Bits(botCellRow_uint_);

			 const __m128i cellsInCols_ = _mm_add_epi8(topCellRow_, _mm_add_epi8(curCellRow_, botCellRow_));

			 currentRemainder_out.cellsCols = (cellsInCols_.m128i_i16[7] & 0b11110000'11110000) >> 4;
			 currentRemainder_out.curCell = (curCellRow_.m128i_i8[15] & 0b11110000) >> 4;

			 const __m128i cellsInCols_sl1_ = shiftLeft1Data(cellsInCols_);
			 const __m128i cellsInCols_sl2_ = shiftLeft1Data(cellsInCols_sl1_);
			 const __m128i cells3by3 =
				 _mm_add_epi8(
					 _mm_add_epi8(
						 _mm_add_epi8(
							 cellsInCols_,
							 cellsInCols_sl1_
						 ),
						 cellsInCols_sl2_
					 ),
					 _mm_srli_si128(
						 _mm_set1_epi16(previousRemainder.cellsCols + (previousRemainder.cellsCols >> 8))
						 , 14
					 )
				 );
			 const __m128i curCellRow = _mm_add_epi8(shiftLeft1Data(curCellRow_), _mm_srli_si128(_mm_set1_epi8(previousRemainder.curCell), 15));

			 const __m128i cellsNeighboursAlive = _mm_sub_epi8(cells3by3, curCellRow);

			 const __m128i aliveCell = _mm_set1_epi8(misc::to_underlying(FieldCell::ALIVE) | (misc::to_underlying(FieldCell::ALIVE) << 4));

			 const __m128i isDead_aliveBits = _mm_andnot_si128(curCellRow, aliveCell); static_assert(
				 (~misc::to_underlying(FieldCell::DEAD)& misc::to_underlying(FieldCell::ALIVE)) == misc::to_underlying(FieldCell::ALIVE) &&
				 (~misc::to_underlying(FieldCell::ALIVE) & misc::to_underlying(FieldCell::ALIVE)) == 0
				 , "need to be able to turn dead into living ones");
			 const __m128i isAlive_aliveBits = curCellRow;

			 const __m128i maskLower = _mm_set1_epi8(0b1111u);
			 const __m128i maskHigher = _mm_set1_epi8(0b1111u << 4);

			 const __m128i cellsNeighboursAliveLower = _mm_and_si128(cellsNeighboursAlive, maskLower);
			 const __m128i is2AliveLower = _mm_and_si128(_mm_cmpeq_epi8(cellsNeighboursAliveLower, _mm_set1_epi8(2u)), maskLower);
			 const __m128i is3AliveLower = _mm_and_si128(_mm_cmpeq_epi8(cellsNeighboursAliveLower, _mm_set1_epi8(3u)), maskLower);

			 const __m128i cellsNeighboursAliveHigher = _mm_and_si128(cellsNeighboursAlive, maskHigher);
			 const __m128i is2AliveHigher = _mm_and_si128(_mm_cmpeq_epi8(cellsNeighboursAliveHigher, _mm_set1_epi8(2u << 4)), maskHigher);
			 const __m128i is3AliveHigher = _mm_and_si128(_mm_cmpeq_epi8(cellsNeighboursAliveHigher, _mm_set1_epi8(3u << 4)), maskHigher);

			 const __m128i is2Alive = _mm_or_si128(is2AliveLower, is2AliveHigher);
			 const __m128i is3Alive = _mm_or_si128(is3AliveLower, is3AliveHigher);

			 const __m128i newGen =
				 _mm_or_si128(
					 _mm_and_si128(
						 isAlive_aliveBits,
						 _mm_or_si128(is2Alive, is3Alive)
					 ),
					 _mm_and_si128(
						 isDead_aliveBits,
						 is3Alive
					 )
				 ); static_assert(misc::to_underlying(FieldCell::DEAD) == 0, "when niether condition is satisfied, resuling bits are 0's");

			 const uint32_t newGen_uint = bits4To32(newGen);
			 return newGen_uint;
	 };

	 const auto isCell = [&grid](int32_t index) -> bool {
		 return grid->isCellAlive(index);
	 };

	 const auto cellAt = [&grid](int32_t index) -> FieldCell {
		 return grid->isCellAlive(index) ? FieldCell::ALIVE : FieldCell::DEAD;
	 };

	 const int32_t startBatch = data->startBatch;
	 const int32_t endBatch = data->endBatch;
	 const int32_t startIndex = startBatch * sizeOfBatch;
	 const int32_t endIndex = endBatch * sizeOfBatch; //may be out of grid bounds. Padding after grid required

	 int32_t i_batch = startBatch;

	 const auto i = [&i_batch]() -> int32_t { return i_batch * sizeOfBatch; };

	 remainder previousRemainder{
		 isCell(i() - 2 - width) +
		 isCell(i() - 2 + 0) +
		 isCell(i() - 2 + width) +
		 (
			 uint16_t(
				 isCell(i() - 1 - width) +
				 isCell(i() - 1 + 0) +
				 isCell(i() - 1 + width))
			 << 8
			 ),
		 isCell(i() - 1 + 0)
	 };

	 for (const uint32_t j_count = 8; (i_batch + j_count) < endBatch;) {
		 for (uint32_t j = 0; j < j_count; ++j, ++i_batch) {
			 const uint32_t newGen = calcNewGen_batch(previousRemainder, i_batch, previousRemainder/*out param*/);

			 grid->setCellsu<Field::FieldPimpl::buffer>(i() - 1, newGen);
		 }

		 if (data->interrupt_flag.load()) return;
	 }

	 for (; i_batch < endBatch; ++i_batch) {
		 const uint32_t newGen = calcNewGen_batch(previousRemainder, i_batch, previousRemainder/*out param*/);

		 grid->setCellsu<Field::FieldPimpl::buffer>(i() - 1, newGen);
	 }

	 if (data->interrupt_flag.load()) return;

	 const uint32_t startRow = startIndex / width;
	 const uint32_t endRow = endIndex / width;

	 {
		 uint8_t //top/cur/bot + first/second/last
			 tf = isCell(-width),
			 ts = isCell(-width + 1),
			 tl = isCell(-width + lastElement),
			 cf = isCell(0),
			 cs = isCell(1),
			 cl = isCell(lastElement);

		 for (uint32_t rowIndex = startRow; rowIndex < endRow; rowIndex++) {
			 const uint32_t row = rowIndex * width;

			 uint8_t
				 bf = isCell(row + width),
				 bs = isCell(row + width + 1),
				 bl = isCell(row + width + lastElement);
			 //first element
			 {
				 const int index = row;
				 const auto curCell = cellAt(index);
				 const unsigned char    topRowNeighbours = tf + ts + tl;
				 const unsigned char bottomRowNeighbours = bf + bs + bl;
				 const unsigned char aliveNeighbours = topRowNeighbours + cl + cs + bottomRowNeighbours;

				 setBufferCellAt(index, fieldCell::nextGeneration(curCell, aliveNeighbours));
			 }

			 tf = cf;
			 ts = cs;
			 tl = cl;
			 cf = bf;
			 cs = bs;
			 cl = bl;
		 }
	 }
	 if (data->interrupt_flag.load()) return;

	 {
		 uint8_t //top/cur/bot + first/last/pre-last
			 tf = isCell(-width),
			 tl = isCell(-width + lastElement),
			 tp = isCell(-width + lastElement - 1),
			 cf = isCell(0),
			 cl = isCell(lastElement),
			 cp = isCell(lastElement - 1);

		 for (uint32_t rowIndex = startRow; rowIndex < endRow; rowIndex++) {
			 const uint32_t row = rowIndex * width;

			 uint8_t
				 bf = isCell(row + width),
				 bp = isCell(row + width + lastElement - 1),
				 bl = isCell(row + width + lastElement);

			 { //last element
				 const int index = row + lastElement;
				 const auto curCell = cellAt(index);
				 const unsigned char    topRowNeighbours = tp + tl + tf;
				 const unsigned char bottomRowNeighbours = bp + bl + bf;
				 const unsigned char aliveNeighbours = topRowNeighbours + cp + cf + bottomRowNeighbours;

				 setBufferCellAt(index, fieldCell::nextGeneration(curCell, aliveNeighbours));
			 }

			 tf = cf;
			 tp = cp;
			 tl = cl;
			 cf = bf;
			 cp = bp;
			 cl = bl;
		 }
	 }
	 if (data->interrupt_flag.load()) return;

	 data->gridUpdate.add(t.elapsedTime());

	 Timer<> t2{};
	 auto& gpuBufferLock_flag = data->gpuBufferLock_flag;
	 bool expected = false;
	 while (!gpuBufferLock_flag.compare_exchange_weak(expected, true)) { expected = false; }// ::std::cout << data->index << ":waiting" << ::std::endl; }

	 data->waiting.add(t2.elapsedTime());

	 const auto bufferP = data->bufferP;
	 glBindBuffer(GL_SHADER_STORAGE_BUFFER, bufferP);
	 glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, bufferP);
	 glBufferSubData(GL_SHADER_STORAGE_BUFFER, data->currentWriteOffset_bytes() + startBatch * sizeOfUint32, (endBatch - startBatch) * sizeOfUint32, &grid->getCellsInt<Field::FieldPimpl::buffer>(startBatch));
	 
	 glUnmapBuffer(GL_SHADER_STORAGE_BUFFER);
	 glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

	 glFinish();
	 data->bufferSend.add(t2.elapsedTime());

	 data->generationUpdated();

	 gpuBufferLock_flag.store(false);
 }



 uint32_t Field::size_bytes() const {
	return misc::roundUpIntTo(size(), sizeOfBatch) / 8;
}

uint32_t* Field::rawData() const {
	return &this->gridPimpl->getCellsInt(0);
}

Field::Field(const uint32_t gridWidth, const uint32_t gridHeight, const size_t numberOfTasks_, const GLuint bufferP_, GLFWwindow* window, bool deployTasks) :
	isFieldGPUBufferWriteOffset{ true }, 
	gridPimpl{ new FieldPimpl(gridWidth, gridHeight) },
	shouldUpdateGrid(true),
	numberOfTasks(numberOfTasks_),
	gridTasks{ new std::unique_ptr<Task<std::unique_ptr<GridData>>>[numberOfTasks_] },
	indecesToBrokenCells{ },
	bufferP{ bufferP_ },
	isStopped{ deployTasks == false }
{
	assert(numberOfTasks >= 1);

	//grid__.median() == 35 ms, bufferSend__.median() == 110 ms.In my case
	//const double fraction = 35.0 / 110.0; // fractionOfWorkTimeToSendData

	const double fraction = 1; // fractionOfWorkTimeToSendData
	static_assert(true, ""/*
	find amountOfWorkT1. given fraction, workload (total batches)
	amountOfWorkT1 — percent of work;

	timeForTask1 = amountOfWorkT1 + amountOfWorkT1 * fraction;
	timeForTask2 = timeForTask1 + timeForTask1 * fraction;
	timeForTask3 = timeForTask2 + timeForTask2 * fraction;
	...


	timeForTask1 + timeForTask2 + timeForTask3 + ... = 1.0;
	amountOfWorkT1 * (1.0 + fraction) + (amountOfWorkT1 * (1.0 + fraction)) * (1.0  + fraction) + ... = workload * (1.0 + fraction);
	amountOfWorkT1 * (1.0 + fraction) + amountOfWorkT1 * (1.0  + fraction)^2 + ...  + amountOfWorkT1 * (1.0  + fraction)^numberOfTasks = workload * (1.0 + fraction);
	amountOfWorkT1 + amountOfWorkT1 * (1.0  + fraction)^1 + ... = workload;
	amountOfWorkT1 * ((1.0 + fraction) + (1.0  + fraction)^2 + ... + (1.0  + fraction)^(numberOfTasks - 1)) = workload;

	amountOfWorkT1 = 1.0 / ((1.0 + fraction) + (1.0  + fraction)^2 + ... + (1.0  + fraction)^(numberOfTasks - 1));

	amountOfWorkT1 = workload / (((1.0 + fraction) * ((1.0 + fraction)^(numberOfTasks-1) - 1)) / fraction)
	amountOfWorkT1 = workload / (((1.0 + fraction) * ((1.0 + fraction)^(numberOfTasks-1) - 1))) * fraction
	*/);

	const auto createGridTask = [window, this](const uint32_t index, const uint32_t startBatch, const uint32_t endBatch) -> void {
		glfwWindowHint(GLFW_VISIBLE, GLFW_FALSE);
		auto* taskWindow = glfwCreateWindow(2, 2, "", NULL, window);
		if (!taskWindow) {
			::std::cout << "task window is null";
			exit(-1);
		}

		gridTasks.get()[index] = std::unique_ptr<Task<std::unique_ptr<GridData>>>(
			new Task<std::unique_ptr<GridData>>{
				threadUpdateGrid,
				std::make_unique<GridData>(
					GridData{
						index,
						this->gridPimpl,
						this->interrupt_flag,
						startBatch,
						endBatch,
						bufferP,
						isFieldGPUBufferWriteOffset,
						size_bytes(),
						gpuBufferLock_flag,
						taskWindow
					}
				)
			}
		);
	};

	uint32_t batchesBefore = 0;
	uint32_t remainingBatches = misc::roundUpIntTo(size(), sizeOfBatch) / sizeOfBatch;

	const double amountOfWorkT1 = remainingBatches / ((1.0 + fraction) * (pow((1.0 + fraction), numberOfTasks - 1) - 1)) * fraction;

	double currentTaskWork = amountOfWorkT1;

	for (uint32_t i = 0; i < numberOfTasks - 1; i++) {
		const uint32_t numberOfBatches = misc::min(remainingBatches, static_cast<uint32_t>(ceil(currentTaskWork)));
		::std::cout << numberOfBatches << std::endl;

		createGridTask(i, batchesBefore, batchesBefore + numberOfBatches);

		batchesBefore    += numberOfBatches;
		remainingBatches -= numberOfBatches;

		currentTaskWork *= (1.0 + fraction);
	}

	{//last
		const uint32_t numberOfBatches = remainingBatches;
		::std::cout << numberOfBatches << std::endl;

		createGridTask(numberOfTasks - 1, batchesBefore, batchesBefore + numberOfBatches);

		batchesBefore    += numberOfBatches;
		remainingBatches -= numberOfBatches;
	}

	deployGridTasks();

	

}
Field::~Field() = default;

uint32_t Field::currentWriteOffset() {
	return size_bytes() * isFieldGPUBufferWriteOffset;
}

uint32_t Field::currentReadOffset() {
	return size_bytes() * !isFieldGPUBufferWriteOffset;
}

void Field::fill(const FieldCell cell) {
	interrupt_flag.store(true);
	waitForGridTasks();
	interrupt_flag.store(false);

	gridPimpl->fill(cell);

	shouldUpdateGrid = true;

	indecesToBrokenCells.clear();
	deployGridTasks();

	//glBindBuffer(GL_SHADER_STORAGE_BUFFER, bufferP);
	//glBufferSubData(GL_SHADER_STORAGE_BUFFER, currentReadOffset() + 0, misc::roundUpIntTo(size(), sizeOfUint32), &gridPimpl->currentCells(0));
	//glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, bufferP);
	//glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
}

FieldCell Field::cellAtIndex(const uint32_t index) const {
	return gridPimpl->cellAt<>(index);
}

void Field::setCellAtIndex(const uint32_t index, FieldCell cell) {
	const auto gridCell{ gridPimpl->cellAt<Field::FieldPimpl::current>(index) };
	if (gridCell != cell) {
		gridPimpl->setCellAt(index, cell);
		//glBindBuffer(GL_SHADER_STORAGE_BUFFER, bufferP);
		//glBufferSubData(GL_SHADER_STORAGE_BUFFER, currentReadOffset() + index, sizeof(), &cell);
		//glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, bufferP);
		//glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
		if (!isStopped) {
			indecesToBrokenCells.push_back(index);

			if (indecesToBrokenCells.size() > (size() / 8)) {
				interrupt_flag.store(true);
				waitForGridTasks();
				interrupt_flag.store(false);
				deployGridTasks();

				indecesToBrokenCells.clear();
			}
		}
		shouldUpdateGrid = true;
	}
}

void Field::updateGeneration() {
	waitForGridTasks();

	if (indecesToBrokenCells.size() > 0) {
		//glBindBuffer(GL_SHADER_STORAGE_BUFFER, bufferP);
		
		for (auto it = indecesToBrokenCells.begin(); it != indecesToBrokenCells.end(); it++) {
			const uint32_t index = *it;

			for (int yo = -1; yo <= 1; yo++) {
				for (int xo = -1; xo <= 1; xo++) {
					int x = ((index + xo) + width()) % width(),
						y = (((index / width()) + yo) + height()) % height();
					int offsetedIndex = x + width() * y;
					const auto newCell = updatedCell(offsetedIndex, gridPimpl);
					gridPimpl->setCellAt<Field::FieldPimpl::buffer>(offsetedIndex, newCell);
					//glBufferSubData(GL_SHADER_STORAGE_BUFFER, currentWriteOffset() + offsetedIndex, sizeof(FieldCell), &newCell);
				}
			}

		}
		indecesToBrokenCells.clear();

		//glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, bufferP);
		//glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
	}

	gridPimpl->prepareNext();
	shouldUpdateGrid = true;
	isFieldGPUBufferWriteOffset = !isFieldGPUBufferWriteOffset;
	interrupt_flag.store(false);

	deployGridTasks();
}

void Field::waitForGridTasks() {
	if (isStopped) return;
	for (uint32_t i = 0; i < numberOfTasks; i++) {
		gridTasks.get()[i]->waitForResult();
	}
}
void Field::deployGridTasks() {
	if (isStopped) return;
	for (uint32_t i = 0; i < numberOfTasks; i++) {
		gridTasks.get()[i]->start();
	}
}

void Field::stopAllGridTasks() {
	interrupt_flag.store(true);
	for (uint32_t i = 0; i < numberOfTasks; i++) {
		gridTasks.get()[i]->waitForResult();
	}
	isStopped = true;
}

void Field::startAllGridTasks() {
	isStopped = false;
	interrupt_flag.store(false);
	for (uint32_t i = 0; i < numberOfTasks; i++) {
		gridTasks.get()[i]->start();
	}
}

uint32_t Field::width() const {
	return gridPimpl->width;
}
uint32_t Field::height() const {
	return gridPimpl->height;
}
uint32_t Field::size() const {
	return gridPimpl->size;
}