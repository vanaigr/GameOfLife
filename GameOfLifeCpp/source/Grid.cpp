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

constexpr size_t sizeOfBatch = sizeof(__m128i); //dependent on CellsGrid fullSize


bool fieldCell::isAlive(const FieldCell cell) {
	return ((misc::to_underlying(cell)) & 0b1);
}
bool fieldCell::isWall(const FieldCell cell) {
	return (misc::to_underlying(cell) >> 4) & 1;
}
bool fieldCell::isDead(const FieldCell cell) {
	return !isWall(cell) && !isAlive(cell);
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
	const uint32_t width;
	const uint32_t height;
	const uint32_t size;
private:
	const uint32_t startPadding;
	const uint32_t endPadding;
	const uint32_t fullSize;
	std::unique_ptr<FieldCell[]> current; 
	std::unique_ptr<FieldCell[]> buffer;

public:
	FieldPimpl(const uint32_t gridWidth, const uint32_t gridHeight) :
		width(gridWidth), height(gridHeight), size(gridWidth* gridHeight),
		startPadding{ misc::roundUpIntTo(gridWidth + 4u, sizeOfBatch) + 15u }, /*
			top row duplicate before first row 
			+ extra bytes for 
				safety and 
				alignment of current row read in 'void threadUpdateGrid(std::unique_ptr<Field::GridData>&)'
		*/
		endPadding{ misc::roundUpIntTo(gridWidth + sizeOfBatch * 4, sizeOfBatch) }, /*
			bottom row duplicate after last row 
			+ extra bytes for 
				memcopy/vector instructions and 
				left/right boundary cells update 
			all in 'void threadUpdateGrid(std::unique_ptr<Field::GridData>&)'
		*/
		fullSize{ startPadding + size + endPadding },
		current{ new FieldCell[fullSize]{ } }, 
		buffer { new FieldCell[fullSize]{ } }
	{};
	~FieldPimpl() = default;

	void prepareNext() {
		current.swap(buffer);
		std::memcpy(&current[0], &(current[size]), sizeof(FieldCell) * startPadding);
		std::memcpy(&current[size + startPadding], &(current[startPadding]), sizeof(FieldCell) * endPadding);
	}

	void fill(const FieldCell cell) {
		std::fill(&current.get()[0], &current.get()[0] + fullSize, cell);
	}

	void fillBuffer(const FieldCell cell) {
		std::fill(&current.get()[0], &current.get()[0] + fullSize, cell);
	}

	FieldCell& cellAt(const int32_t index) const {
		return current.get()[index + startPadding];
	}

	FieldCell& bufferCellAt(const int32_t index) const {
		return buffer.get()[index + startPadding];
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

	 uint32_t currentWriteOffset() {
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
	const FieldCell* grid = &cellsGrid->cellAt(0);
	const auto width = cellsGrid->width;
	const auto height = cellsGrid->height;
	FieldCell cell = grid[index];

	if (cell != FieldCell::WALL) {
		uint32_t aliveNeighbours = 0;

		for (int yo = -1; yo <= 1; yo++) {
			for (int xo = -1; xo <= 1; xo++) {
				if (xo != 0 || yo != 0) {
					int x = ((index + xo) + width) % width,
						y = (((index / width) + yo) + height) % height;
					int offsetedIndex = x + width * y;
					if (fieldCell::isAlive(cellsGrid->cellAt(offsetedIndex)))
						aliveNeighbours++;
				}
			}
		}

		return fieldCell::nextGeneration(cell, aliveNeighbours);
	}

	return cell;
}


static void threadUpdateGrid(std::unique_ptr<Field::GridData>& data) {
	Timer<> t{};
	auto& grid = data->grid;
	const int32_t width = static_cast<int32_t>(grid->width);
	const int32_t height = static_cast<int32_t>(grid->height);
	const int lastElement = width - 1;

	GLFWwindow* window = data->offscreen_context;
	if(window != glfwGetCurrentContext())glfwMakeContextCurrent(window);

	const auto setBufferCellAt = [&grid](uint32_t index, FieldCell cell) -> void { grid->bufferCellAt(index) = cell; }; 
	const auto cellAt = [&grid](int32_t index) -> FieldCell& { return grid->cellAt(index); };
	const auto isCell = [&cellAt](int32_t index) -> bool { return fieldCell::isAlive(cellAt(index)); };

	struct remainder {
		uint16_t cellsCols;
		uint8_t curCell;
	};

	const auto calcNewGenBatch128 = [&grid](
		const remainder previousRemainder,
		uint32_t index,
		remainder& currentRemainder_out
		) -> __m128i {
		static_assert(sizeOfBatch > 1, "sizeOfBatch = 0, 1 doesnt make sense");
		const auto cellAt = [&grid](uint32_t index) -> FieldCell& { return grid->cellAt(index); };
		const auto width = grid->width;

		static_assert(sizeof FieldCell == 1, "size of field cell must be 1 byte");
		const __m128i topCellRow_ = _mm_loadu_si128(reinterpret_cast<__m128i*>(&cellAt(index - width)));
		const __m128i curCellRow_ = _mm_load_si128 (reinterpret_cast<__m128i*>(&cellAt(index        )));
		const __m128i botCellRow_ = _mm_loadu_si128(reinterpret_cast<__m128i*>(&cellAt(index + width)));

		//if ((uintptr_t(&cellAt(index)) & 0xf) != 0) std::cout << "(sick!)" << (uintptr_t(&cellAt(index)) & 0xf) << std::endl;

		const __m128i cellsInCols_ = _mm_add_epi8(topCellRow_, _mm_add_epi8(curCellRow_, botCellRow_));

		currentRemainder_out.cellsCols = cellsInCols_.m128i_i16[7];
		currentRemainder_out.curCell   = curCellRow_ .m128i_i8 [15];


		const __m128i cells3by3 =
			_mm_add_epi8(
				_mm_add_epi8(
					_mm_slli_si128(cellsInCols_, 1),
					_mm_add_epi8(
						cellsInCols_,
						_mm_slli_si128(cellsInCols_, 2)
					)
				),
				_mm_srli_si128(
					_mm_set1_epi16(previousRemainder.cellsCols + (previousRemainder.cellsCols >> 8))
					, 14
				)
			); static_assert(
				static_cast<uint16_t>(misc::to_underlying(FieldCell::WALL)) * 9 == misc::to_underlying(FieldCell::WALL) * 9 &&
				static_cast<uint16_t>(misc::to_underlying(FieldCell::ALIVE)) * 9 == misc::to_underlying(FieldCell::ALIVE) * 9 &&
				(((misc::to_underlying(FieldCell::WALL)) * 9) & (misc::to_underlying(FieldCell::ALIVE) * 9)) == 0
				, "when adding 9 cells together there should be no overflows, and no inersections between number of alive cells and number of walls"

				);
		const __m128i curCellRow = _mm_add_epi8(_mm_slli_si128(curCellRow_, 1), _mm_srli_si128(_mm_set1_epi8(previousRemainder.curCell), 15));
		const __m128i cellsInCols = _mm_add_epi8(_mm_slli_si128(cellsInCols_, 1), _mm_srli_si128(_mm_set1_epi8(previousRemainder.cellsCols >> 8), 15));


		/*std::cout << "here:" << (_mm_movemask_epi8(
			_mm_cmpeq_epi8(
				_mm_srli_si128( // >> 14*8
					_mm_set1_epi16(previousRemainder.cellsCols + (previousRemainder.cellsCols >> 8)) // most significant -> [ ...-higher_byte-lower_byte ] <-least significant
					, 14
				) // most significant -> [ ...0-0-higher_byte-lower_byte ] <-least significant
				,_mm_setr_epi16(previousRemainder.cellsCols + (previousRemainder.cellsCols >> 8), 0, 0, 0, 0, 0, 0, 0)// most significant -> [ higher_byte-lower_byte-0-0 ... 0-0-0-0 ] <-least significant
			)
		) ==
			_mm_movemask_epi8(
				_mm_set1_epi32(0xffffffff)
			)) << std::endl;*/
			//printed "here:1"
			//why setr and not just set?
			/* actualy:
			intel intrinsics website: __m128i _mm_setr_epi16(short e7 , short e6 , short e5 , short e4 , short e3 , short e2 , short e1 , short e0 )
			visual studio hint      : __m128i _mm_setr_epi16(short _W0, short _W1, short _W2, short _W3, short _W4, short _W5, short _W6, short _W7)

			action performed:
			dst[15:0] := e7
			dst[31:16] := e6
			dst[47:32] := e5
			dst[63:48] := e4
			dst[79:64] := e3
			dst[95:80] := e2
			dst[111:96] := e1
			dst[127:112] := e0

			why intel's "e7" corresponds to hint's "_W0"?
			*/


		static_assert(((
			(misc::to_underlying(FieldCell::ALIVE) + misc::to_underlying(FieldCell::WALL)) * 8)
			& 0b0000'1111) == misc::to_underlying(FieldCell::ALIVE) * 8, "we need to be able to remove all the walls from 3by3 sum of cells");
		const __m128i lowerFourMask = _mm_set1_epi8(0b1111);
		const __m128i cellsNeighboursAlive = _mm_and_si128(_mm_sub_epi8(cells3by3, curCellRow), lowerFourMask);
		
		const __m128i wallCell = _mm_set1_epi8(misc::to_underlying(FieldCell::WALL));
		const __m128i deadCell = _mm_set1_epi8(misc::to_underlying(FieldCell::DEAD));
		const __m128i aliveCell = _mm_set1_epi8(misc::to_underlying(FieldCell::ALIVE));
		static_assert((misc::to_underlying(FieldCell::ALIVE)& misc::to_underlying(FieldCell::WALL)) == 0, "must not have intersections");


		const __m128i isWall = _mm_cmpeq_epi8(curCellRow, wallCell);
		const __m128i wallsBit = _mm_and_si128(curCellRow, wallCell);
		const __m128i isDead_lowerBit = _mm_andnot_si128(curCellRow, aliveCell); static_assert(
			(~misc::to_underlying(FieldCell::DEAD)& misc::to_underlying(FieldCell::ALIVE)) == misc::to_underlying(FieldCell::ALIVE) &&
			(~misc::to_underlying(FieldCell::ALIVE) & misc::to_underlying(FieldCell::ALIVE)) == 0
			, "need to be able to turn dead into living ones");
		const __m128i isAlive_lowerBit = curCellRow;
		const __m128i is2Alive = _mm_cmpeq_epi8(cellsNeighboursAlive, _mm_set1_epi8(2));
		const __m128i is3Alive = _mm_cmpeq_epi8(cellsNeighboursAlive, _mm_set1_epi8(3));

		const __m128i newGen =
			_mm_or_si128(
				wallsBit,
				_mm_andnot_si128(
					isWall,
					_mm_or_si128(
						_mm_and_si128(
							isAlive_lowerBit,
							_mm_or_si128(is2Alive, is3Alive)
						),
						_mm_and_si128(
							isDead_lowerBit,
							is3Alive
						)
					)
				)
			);

		return newGen;
	};


	const int32_t startIndex = data->startBatch * sizeOfBatch;
	const int32_t endIndex   = data->endBatch * sizeOfBatch; //may be out of grid bounds. Padding after grid required

	int32_t i = startIndex;

	static_assert(sizeOfBatch > 2, "");
	remainder previousRemainder{
		misc::to_underlying(cellAt(i - 2 - width)) +
		misc::to_underlying(cellAt(i - 2 + 0)) +
		misc::to_underlying(cellAt(i - 2 + width)) +
		(
			uint16_t( misc::to_underlying(cellAt(i - 1 - width)) +
			  misc::to_underlying(cellAt(i - 1 + 0)) +
			  misc::to_underlying(cellAt(i - 1 + width)) ) 
			<< 8
		),
		misc::to_underlying(cellAt(i - 1 + 0))
	};

	for (const uint32_t j_count = 8; (i + sizeOfBatch * j_count) < endIndex;) {
		for (uint32_t j = 0; j < j_count; j++, i += sizeOfBatch) {
			const __m128i newGen = calcNewGenBatch128(previousRemainder, i + 1, previousRemainder/*out param*/);

			std::memcpy(&grid->bufferCellAt(i), ((FieldCell*)(&newGen)), sizeof(FieldCell) * sizeOfBatch);
		}

		if (data->interrupt_flag.load()) return;
	}

	for (; i < endIndex; i += sizeOfBatch) {
		const __m128i newGen = calcNewGenBatch128(previousRemainder, i + 1, previousRemainder/*out param*/);

		std::memcpy(&grid->bufferCellAt(i), ((FieldCell*)(&newGen)), sizeof(FieldCell) * sizeOfBatch);
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
	glBufferSubData(GL_SHADER_STORAGE_BUFFER, data->currentWriteOffset() + startIndex, sizeof(FieldCell) * (misc::min<uint32_t>(endIndex, grid->size) - startIndex), &grid->bufferCellAt(startIndex));

	//int32_t packIndex = startIndex;
	//static_assert(
	//	misc::to_underlying(FieldCell::WALL) == 0b0001'0000 && 
	//	misc::to_underlying(FieldCell::ALIVE) == 0b0000'0001 && 
	//	misc::to_underlying(FieldCell::DEAD) == 0,
	//	"required values for this kind of packing");
	//for (; packIndex < endIndex; packIndex += sizeOfBatch * 4 /*safe because of 16u*4 padding at the end of the grid*/) {
	//	__m128i packedCells = _mm_set1_epi8(0);
	//	for (int32_t j = 0; j < 4; j++)
	//		packedCells =
	//			_mm_or_si128(
	//				packedCells,
	//				_mm_slli_epi16(//because there is no _mm_slli_epi8
	//					_mm_load_si128(reinterpret_cast<__m128i*>(&cellAt(packIndex + j * sizeOfBatch)))
	//					, j
	//				) 
	//		);

	//	glBufferSubData(GL_SHADER_STORAGE_BUFFER, data->currentOffset() + startIndex, sizeOfBatch, &packedCells);
	//}
	glUnmapBuffer(GL_SHADER_STORAGE_BUFFER);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

	glFinish();
	data->bufferSend.add(t2.elapsedTime());

	data->generationUpdated();

	gpuBufferLock_flag.store(false);
}


FieldCell* Field::grid() {
	shouldUpdateGrid = false;
	return &this->gridPimpl->cellAt(0);
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
						size(),
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
	return size() * isFieldGPUBufferWriteOffset;
}

uint32_t Field::currentReadOffset() {
	return size() * !isFieldGPUBufferWriteOffset;
}

void Field::fill(const FieldCell cell) {
	interrupt_flag.store(true);
	waitForGridTasks();
	interrupt_flag.store(false);

	gridPimpl->fill(cell);

	shouldUpdateGrid = true;

	indecesToBrokenCells.clear();
	deployGridTasks();

	glBindBuffer(GL_SHADER_STORAGE_BUFFER, bufferP);
	glBufferSubData(GL_SHADER_STORAGE_BUFFER, currentReadOffset() + 0, sizeof(FieldCell) * size(), &gridPimpl->cellAt(0));
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, bufferP);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
}

FieldCell Field::cellAtIndex(const uint32_t index) const {
	return gridPimpl->cellAt(index);
}

void Field::setCellAtIndex(const uint32_t index, FieldCell cell) {
	auto& gridCell = gridPimpl->cellAt(index);
	if (gridCell != cell) {
		gridCell = cell;
		glBindBuffer(GL_SHADER_STORAGE_BUFFER, bufferP);
		glBufferSubData(GL_SHADER_STORAGE_BUFFER, currentReadOffset() + index, sizeof(FieldCell), &cell);
		glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, bufferP);
		glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
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
		glBindBuffer(GL_SHADER_STORAGE_BUFFER, bufferP);
		
		for (auto it = indecesToBrokenCells.begin(); it != indecesToBrokenCells.end(); it++) {
			const uint32_t index = *it;

			for (int yo = -1; yo <= 1; yo++) {
				for (int xo = -1; xo <= 1; xo++) {
					int x = ((index + xo) + width()) % width(),
						y = (((index / width()) + yo) + height()) % height();
					int offsetedIndex = x + width() * y;
					const auto newCell = updatedCell(offsetedIndex, gridPimpl);
					gridPimpl->bufferCellAt(offsetedIndex) = newCell;
					glBufferSubData(GL_SHADER_STORAGE_BUFFER, currentWriteOffset() + offsetedIndex, sizeof(FieldCell), &newCell);
				}
			}

		}
		indecesToBrokenCells.clear();

		glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, bufferP);
		glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
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