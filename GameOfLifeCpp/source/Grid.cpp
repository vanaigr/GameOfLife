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

static const size_t sizeOfUint32 = sizeof uint32_t;
constexpr uint32_t batchSize = 32;// sizeof(__m128i); //dependent on CellsGrid fullSize

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


 class Field::FieldPimpl {
 public:
	 struct current {};
	 struct buffer  {};
public:
	const uint32_t width_grid;
	const uint32_t width_actual;
	const uint32_t width_int;
	const uint32_t height;
	const uint32_t size_grid;
	const uint32_t size_actual;
	const uint32_t size_int;
private:
	const uint32_t startPadding_int;
	const uint32_t endPadding_int;
	const uint32_t fullSize_int;
	std::unique_ptr<uint32_t[]> current_; 
	std::unique_ptr<uint32_t[]> buffer_;

public:
	FieldPimpl(const uint32_t gridWidth, const uint32_t gridHeight) :
		width_grid{ gridWidth },
		width_actual{ misc::roundUpIntTo(width_grid, batchSize) },
		width_int{ misc::intDivCeil(width_actual, batchSize) },

		height{ gridHeight },

		size_grid(width_grid * height),
		size_actual{ width_actual * height },
		size_int{ misc::intDivCeil(width_actual * height, batchSize) },

		startPadding_int{ width_int }, /*
			top row duplicate before the first row
		*/
		endPadding_int{ misc::intDivCeil(width_grid, batchSize) + 1 }, /*
			bottom row duplicate after the last row
			+ extra bytes for
				memcopy/vector instructions and
				left/right boundary cells update
			all in 'void threadUpdateGrid(std::unique_ptr<Field::GridData>&)'
		*/
		fullSize_int{ startPadding_int + size_int + endPadding_int },
		current_{ new uint32_t[fullSize_int]{ } },
		buffer_{ new uint32_t[fullSize_int]{ } }
	{ };
	~FieldPimpl() = default;

	void prepareNext() {
		current_.swap(buffer_);
		std::memcpy(&current_[0], &(current_[size_int]), width_int * sizeOfUint32); //start pading
		std::memcpy(&current_[startPadding_int + size_int], &(current_[startPadding_int]), width_int * sizeOfUint32); //end pading
		
	}

	template<typename gridType = current>
	void fill(const FieldCell cell) {
		auto* grid{ getGrid<gridType>() };
		const auto val = ~0u * fieldCell::isAlive(cell);
		std::fill(&grid[0], &grid[0] + fullSize_int, val);
	}

	template<typename gridType = current>
	bool isCellAlive(const int32_t index) const {
		auto* grid{ getGrid<gridType>() };
		const auto row = index / int32_t(width_grid);
		const auto col = misc::mod(index, width_grid);

		const auto col_int = col / batchSize;
		const auto shift = col % batchSize;

		return (grid[startPadding_int + row * width_int + col_int] >> shift) & 0b1;
	}

	template<typename gridType = current>
	bool isCellAlive_actual(const int32_t index) const {
		auto* grid{ getGrid<gridType>() };

		const auto offset = index / int32_t(batchSize);
		const auto shift = misc::mod(index, batchSize);

		return (grid[startPadding_int + offset] >> shift) & 0b1;
	}

	template<typename gridType = current>
	FieldCell cellAt(const int32_t index) const {
		return isCellAlive<gridType>(index) ? FieldCell::ALIVE : FieldCell::DEAD;
	}

	template<typename gridType = current>
	void setCellAt(const int32_t index, const FieldCell cell) const {
		auto* grid{ getGrid<gridType>() };

		const auto row = index / width_grid;
		const auto col = misc::mod(index, width_grid);

		const auto col_int = col / batchSize;
		const auto shift = misc::mod(col, batchSize);

		auto& cur{ grid[startPadding_int + row * width_int + col_int] };

		cur = (cur & ~(0b1 << shift)) | (fieldCell::isAlive(cell) << shift);
	}/* tested */

	template<typename gridType = current>
	uint32_t& getCellsActual_int(int32_t index_actual_int) {
		return getGrid<gridType>()[index_actual_int + startPadding_int];
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
	 const auto width_grid = cellsGrid->width_grid;
	 const auto height = cellsGrid->height;
	 bool cell = cellsGrid->isCellAlive(index);

	 uint32_t aliveNeighbours = 0;

	 for (int yo = -1; yo <= 1; yo++) {
		 for (int xo = -1; xo <= 1; xo++) {
			 if (xo != 0 || yo != 0) {
				 int x = ((index + xo) + width_grid) % width_grid,
					 y = (((index / width_grid) + yo) + height) % height;
				 int offsetedIndex = x + width_grid * y;
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
	 const int32_t width_grid = static_cast<int32_t>(grid->width_grid);
	 const int32_t width_actual = static_cast<int32_t>(grid->width_actual);
	 const int32_t height = static_cast<int32_t>(grid->height);
	 const int lastElement = width_grid - 1;

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

			 static_assert(batchSize > 1, "sizeOfBatch = 0, 1 doesnt make sense");
			 const auto width_grid = grid->width_grid;
			 const auto width_actual = grid->width_actual;
			 const auto width_int = grid->width_int;

			 const uint32_t topCellRow_uint_{ grid->getCellsActual_int(index_batch - width_int) };
			 const uint32_t curCellRow_uint_{ grid->getCellsActual_int(index_batch) };
			 const uint32_t botCellRow_uint_{ grid->getCellsActual_int(index_batch + width_int) };

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

			 const __m128i aliveCell = _mm_set1_epi8(0b00010001);

			 const __m128i isDead_aliveBits = _mm_andnot_si128(curCellRow, aliveCell); 
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
				 );

			 const uint32_t newGen_uint = bits4To32(newGen);
			 return newGen_uint;
	 };

	 const auto isCell = [&grid](int32_t index) -> bool {
		 return grid->isCellAlive(index);
	 };

	 const auto isCell_actual = [&grid](int32_t index) -> bool {
		 return grid->isCellAlive_actual(index);
	 };

	 const auto cellAt = [&grid](int32_t index) -> FieldCell {
		 return grid->isCellAlive(index) ? FieldCell::ALIVE : FieldCell::DEAD;
	 };

	 const int32_t startBatch = data->startBatch;
	 const int32_t endBatch = data->endBatch;
	 const int32_t startIndex = startBatch * batchSize;
	 const int32_t endIndex = endBatch * batchSize; //may be out of grid bounds. Padding after grid required

	 int32_t i_batch = startBatch;

	 const auto i_a = [=]() -> int32_t { 
		 return i_batch * batchSize; 
	 };

	 remainder previousRemainder{
		 isCell_actual(i_a() - 1 - width_actual) +
		 isCell_actual(i_a() - 1 + 0) +
		 isCell_actual(i_a() - 1 + width_actual) +
		 (
			 uint16_t(
				 isCell_actual(i_a() - width_actual) +
				 isCell_actual(i_a() + 0) +
				 isCell_actual(i_a() + width_actual))
			 << 8
			 ),
		 isCell(i_a())
	 };

	 uint64_t newGenWindow = 0;
	 for (const uint32_t j_count = 8; (i_batch + j_count) < endBatch + 1;) { //endIndex + 1 is because of last celll is not updated
		 for (uint32_t j = 0; j < j_count; ++j, ++i_batch) {
			 const uint32_t newGen = calcNewGen_batch(previousRemainder, i_batch, previousRemainder/*out param*/);

			 newGenWindow = (newGenWindow >> 32) | (uint64_t(newGen) << 31);
			 grid->getCellsActual_int<Field::FieldPimpl::buffer>(i_batch - 1) = uint32_t(newGenWindow);
		 }

		 if (data->interrupt_flag.load()) return;
	 }

	 for (; i_batch < endBatch + 1; ++i_batch) {
		 const uint32_t newGen = calcNewGen_batch(previousRemainder, i_batch, previousRemainder/*out param*/);

		 newGenWindow = (newGenWindow >> 32) | (uint64_t(newGen) << 31);
		 grid->getCellsActual_int<Field::FieldPimpl::buffer>(i_batch - 1) = uint32_t(newGenWindow);
	 }


	 if (data->interrupt_flag.load()) return;

	 const uint32_t startRow = startIndex / width_actual;
	 const uint32_t endRow = misc::min<uint32_t>(misc::intDivCeil(endIndex + 1, width_actual), height);
	 {
		 uint8_t //top/cur/bot + first/second/last
			 tf = isCell(-width_grid),
			 ts = isCell(-width_grid + 1),
			 tl = isCell(-width_grid + lastElement),
			 cf = isCell(0),
			 cs = isCell(1),
			 cl = isCell(lastElement);

		 for (uint32_t rowIndex = startRow; rowIndex < endRow; rowIndex++) {
			 const uint32_t row = rowIndex * width_grid;

			 uint8_t
				 bf = isCell(row + width_grid),
				 bs = isCell(row + width_grid + 1),
				 bl = isCell(row + width_grid + lastElement);
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
			 tf = isCell(-width_grid),
			 tl = isCell(-width_grid + lastElement),
			 tp = isCell(-width_grid + lastElement - 1),
			 cf = isCell(0),
			 cl = isCell(lastElement),
			 cp = isCell(lastElement - 1);

		 for (uint32_t rowIndex = startRow; rowIndex < endRow; rowIndex++) {
			 const uint32_t row = rowIndex * width_grid;

			 uint8_t
				 bf = isCell(row + width_grid),
				 bp = isCell(row + width_grid + lastElement - 1),
				 bl = isCell(row + width_grid + lastElement);

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
	 glBufferSubData(GL_SHADER_STORAGE_BUFFER, data->currentWriteOffset_bytes() + startBatch * sizeOfUint32, (endBatch - startBatch) * sizeOfUint32, &grid->getCellsActual_int<Field::FieldPimpl::buffer>(startBatch));
	 
	 glUnmapBuffer(GL_SHADER_STORAGE_BUFFER);
	 glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

	 glFinish();
	 data->bufferSend.add(t2.elapsedTime());

	 data->generationUpdated();

	 gpuBufferLock_flag.store(false);
 }

 uint32_t Field::width_actual() const {
	 return gridPimpl->width_actual;
 }

 uint32_t Field::size_bytes() const {
	return gridPimpl->size_int * sizeOfUint32;
}

 uint32_t Field::size_actual() const {
	 return gridPimpl->size_actual;
 }

uint32_t* Field::rawData() const {
	return &this->gridPimpl->getCellsActual_int(0);
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

	const double fraction = .1; // fractionOfWorkTimeToSendData
	static_assert(true, ""/*
	find amountOfWorkT1. given fraction, workload (total batches)
	amountOfWorkT1 — percent of work;

	timeForTask1 = amountOfWorkT1 + amountOfWorkT1 * fraction;
	timeForTask2 = timeForTask1 + timeForTask1 * fraction;
	timeForTask3 = timeForTask2 + timeForTask2 * fraction;
	...

	                                                            v  total time  
	timeForTask1 + timeForTask2 + timeForTask3 + ... = workload * (1.0 + fraction);
	amountOfWorkT1 * (1.0 + fraction) 
		+ (amountOfWorkT1 * (1.0 + fraction)) * (1.0  + fraction) 
		+ ...                                                       = workload * (1.0 + fraction);
	amountOfWorkT1 * (1.0 + fraction) + amountOfWorkT1 * (1.0  + fraction)^2 + ...  + amountOfWorkT1 * (1.0  + fraction)^numberOfTasks = workload * (1.0 + fraction);

	amountOfWorkT1 * ((1.0 + fraction) + (1.0  + fraction)^2 + ... + (1.0  + fraction)^numberOfTasks) = workload * (1.0 + fraction);

	amountOfWorkT1 = workload / ( ( (1.0 + fraction)^numberOfTasks - 1 ) / fraction );
	amountOfWorkT1 = workload /   ( (1.0 + fraction)^numberOfTasks - 1 ) * fraction  ;
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
	uint32_t remainingBatches = misc::intDivCeil(size_actual(), batchSize);

	const double amountOfWorkT1 = remainingBatches / (pow((1.0 + fraction), numberOfTasks) - 1) * fraction;

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
	return gridPimpl->width_grid;
}
uint32_t Field::height() const {
	return gridPimpl->height;
}
uint32_t Field::size() const {
	return gridPimpl->size_grid;
}