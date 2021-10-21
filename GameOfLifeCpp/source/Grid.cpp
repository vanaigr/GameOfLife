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
	const uint32_t padding;
	const uint32_t fullSize;
	std::unique_ptr<FieldCell[]> current; 
	std::unique_ptr<FieldCell[]> buffer;

public:
	FieldPimpl(const uint32_t gridWidth, const uint32_t gridHeight) :
		width(gridWidth), height(gridHeight), size(gridWidth* gridHeight),
		padding{ gridWidth },
		fullSize{ padding + size + padding + 16u },//+2 duplicate rows (top row duplicate before first row, bottom row duplicate after last row) + extra bytes for memcopy/vector instruction in thread grid update
		current{ new FieldCell[fullSize]{ FieldCell::DEAD } }, 
		buffer{ new FieldCell[fullSize]{ FieldCell::DEAD } }
	{};
	~FieldPimpl() = default;

	void prepareNext() {
		current.swap(buffer);
		std::memcpy(&current[0], &(current[size]), sizeof(FieldCell) * padding);
		std::memcpy(&current[size + padding], &(current[padding]), sizeof(FieldCell) * padding);
	}

	void fill(const FieldCell cell) {
		std::fill(&current.get()[0], &current.get()[0] + fullSize, cell);
	}

	void fillBuffer(const FieldCell cell) {
		std::fill(&current.get()[0], &current.get()[0] + fullSize, cell);
	}

	FieldCell& cellAt(const int32_t index) const {
		return current.get()[index + padding];
	}

	FieldCell& bufferCellAt(const int32_t index) const {
		return buffer.get()[index + padding];
	}
};

struct Field::GridData {
	MedianCounter gridUpdate{}, waiting{}, bufferSend{};
	unsigned int grid__iteration = 0;

	uint32_t index;

	std::unique_ptr<FieldPimpl>& grid;
	std::atomic_bool& interrupt_flag;
	std::atomic_bool& gpuBufferLock_flag;
	uint32_t startRow;
	uint32_t rowCount;

	GLFWwindow* offscreen_context;
	GLuint bufferP;
	bool& isOffset;
	uint32_t offset;

	uint32_t currentOffset() {
		return offset * isOffset;
	}

public:
	GridData(
		uint32_t index_,
		std::unique_ptr<FieldPimpl>& grid_,
		std::atomic_bool& interrupt_flag_,
		uint32_t startRow_,
		uint32_t rowCount_,
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
		startRow(startRow_),
		rowCount(rowCount_),
		offscreen_context(offscreen_context_),
		bufferP(bufferP_),
		isOffset(isOffset_),
		offset(offset_)
	{}
};

//struct Field::PackedGridData {
//	std::unique_ptr<FieldPimpl>& grid;
//	std::shared_ptr<uint32_t[]>& packedGrid;
//
//	uint32_t startCell;
//	uint32_t cellsCount;
//};

FieldCell updatedCell(const int32_t index, const std::unique_ptr<Field::FieldPimpl>& cellsGrid) {
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

void threadUpdateGrid(std::unique_ptr<Field::GridData>& data) {
	Timer<> t{};
	auto& grid = data->grid;
	const int32_t width = static_cast<int32_t>(grid->width);
	const int32_t height = static_cast<int32_t>(grid->height);
	const int lastElement = width - 1;
	const auto rowCount = data->rowCount;
	const auto startRow = data->startRow;

	GLFWwindow* window = data->offscreen_context;
	glfwMakeContextCurrent(window);

	const auto setBufferCellAt = [&grid](uint32_t index, FieldCell cell) -> void { grid->bufferCellAt(index) = cell; }; // { return data.grid[data.gridO.normalizeIndex(index + (data.startRow * data.gridO.width))]; };
	const auto cellAt = [&grid](int32_t index) -> FieldCell& { return grid->cellAt(index); }; // { return data.grid[data.gridO.normalizeIndex(index + (data.startRow * data.gridO.width))]; };
	const auto isCell = [&cellAt](int32_t index) -> bool { return fieldCell::isAlive(cellAt(index)); };

	const auto calcNewGenBatch128 = [&grid](
		uint32_t index
		) -> __m128i {
		const auto cellAt = [&grid](uint32_t index) -> FieldCell& { return grid->cellAt(index); }; // { return data.grid[data.gridO.normalizeIndex(index + (data.startRow * data.gridO.width))]; };
		const auto width = grid->width;

		const __m128i topCellRow = _mm_loadu_si128(reinterpret_cast<__m128i*>(&cellAt(index - width)));
		const __m128i curCellRow = _mm_loadu_si128(reinterpret_cast<__m128i*>(&cellAt(index)));
		const __m128i botCellRow = _mm_loadu_si128(reinterpret_cast<__m128i*>(&cellAt(index + width)));

		const __m128i lowerFourMask = _mm_set1_epi8(0b1111);
		const __m128i wallCell = _mm_set1_epi8(misc::to_underlying(FieldCell::WALL));
		const __m128i deadCell = _mm_set1_epi8(misc::to_underlying(FieldCell::DEAD));
		const __m128i aliveCell = _mm_set1_epi8(misc::to_underlying(FieldCell::ALIVE));

		const __m128i cellsInRows = _mm_add_epi8(topCellRow, _mm_add_epi8(curCellRow, botCellRow));
		const __m128i cells3by3Neighbours =
			_mm_add_epi8(
				cellsInRows,
				_mm_add_epi8(
					_mm_slli_si128(cellsInRows, 1), // << 8 
					_mm_srli_si128(cellsInRows, 1)  // >> 8
				)
			);

		const __m128i cellsActualNeighboursAlive = _mm_and_si128(_mm_sub_epi8(cells3by3Neighbours, curCellRow), lowerFourMask);


		const __m128i walls    = _mm_cmpeq_epi8(curCellRow, wallCell);
		const __m128i dead     = _mm_cmpeq_epi8(curCellRow, deadCell);
		const __m128i alive    = _mm_cmpeq_epi8(curCellRow, aliveCell);
		const __m128i is2Alive = _mm_cmpeq_epi8(cellsActualNeighboursAlive, _mm_set1_epi8(2));
		const __m128i is3Alive = _mm_cmpeq_epi8(cellsActualNeighboursAlive, _mm_set1_epi8(3));

		const __m128i newGen =
			_mm_or_si128(
				_mm_and_si128(walls, wallCell),
				_mm_and_si128(_mm_or_si128(
					_mm_and_si128(
						alive,
						_mm_or_si128(is2Alive, is3Alive)
					),
					_mm_and_si128(
						dead,
						is3Alive
					)
				), aliveCell)
			);

		return newGen;
	};

	const uint32_t startIndex = startRow * width;
	const uint32_t endIndex = (startRow + rowCount) * width;

	const size_t sizeOfBatch = sizeof(__m128i); //dependent on CellsGrid fullSize
	static_assert(sizeOfBatch > 2, "sizeOfBatch of 0, 1, 2 doesnt make sense here");

	/*const uint32_t horisontalBatches = (width + ((sizeOfBatch - 2) - 1)) / (sizeOfBatch - 2);

	for (uint32_t batch = 0; batch < horisontalBatches - 1; batch++) {
		__m128i topCellRow = _mm_loadu_si128(reinterpret_cast<__m128i*>(&cellAt(0 - width)));
		__m128i curCellRow = _mm_loadu_si128(reinterpret_cast<__m128i*>(&cellAt(0)));
		for (uint32_t row = startRow; row < startRow + rowCount; row++) {
			uint32_t index = row * width + batch * (sizeOfBatch - 2);
			const auto tuple = calcNewGenBatch128(index, topCellRow, curCellRow);
			const __m128i newGen = std::get<0>(tuple);
			topCellRow = std::get<1>(tuple);
			curCellRow = std::get<2>(tuple);

			std::memcpy(&grid->bufferCellAt(index + 1), ((FieldCell*)(&newGen)) + 1, sizeof(FieldCell) * (sizeOfBatch - 2));
		}
		if (data->interrupt_flag.load()) return;
	}

	__m128i topCellRow = _mm_loadu_si128(reinterpret_cast<__m128i*>(&cellAt(0 - width)));
	__m128i curCellRow = _mm_loadu_si128(reinterpret_cast<__m128i*>(&cellAt(0)));
	{ //last row of batches 
		for (uint32_t row = startRow; row < startRow + rowCount - 1; row++) {
			uint32_t index = row * width + (horisontalBatches - 1) * (sizeOfBatch - 2);
			//const __m128i newGen = calcNewGenBatch128(index);
			const auto tuple = calcNewGenBatch128(index, topCellRow, curCellRow);
			const __m128i newGen = std::get<0>(tuple);
			topCellRow = std::get<1>(tuple);
			curCellRow = std::get<2>(tuple);

			std::memcpy(&grid->bufferCellAt(index + 1), ((FieldCell*)(&newGen)) + 1, sizeof(FieldCell) * (sizeOfBatch - 2));
		}
		if (data->interrupt_flag.load()) return;
	}

	{ //last batch
		uint32_t index = (startRow + rowCount - 1) * width + (horisontalBatches - 1) * (sizeOfBatch - 2);
		//const __m128i newGen = calcNewGenBatch128(index);
		const auto tuple = calcNewGenBatch128(index, topCellRow, curCellRow);
		const __m128i newGen = std::get<0>(tuple);
		topCellRow = std::get<1>(tuple);
		curCellRow = std::get<2>(tuple);

		const uint32_t remainingCells = width - (horisontalBatches-1) * (sizeOfBatch - 2);

		std::memcpy(&grid->bufferCellAt(index + 1), ((FieldCell*)(&newGen)) + 1, sizeof(FieldCell) * remainingCells);
	}*/

	uint32_t i = startIndex;


	for (; (i + sizeOfBatch - 2) < endIndex; i += sizeOfBatch - 2) {
		const __m128i newGen = calcNewGenBatch128(i);

		std::memcpy(&grid->bufferCellAt(i + 1), ((FieldCell*)(&newGen)) + 1, sizeof(FieldCell) * (sizeOfBatch - 2));

		if (data->interrupt_flag.load()) return;
	}

	if (i != endIndex) {
		uint32_t remainingCells = endIndex - i;

		const __m128i newGen = calcNewGenBatch128(i);

		std::memcpy(&grid->bufferCellAt(i + 1), ((FieldCell*)(&newGen)) + 1, sizeof(FieldCell) * (remainingCells - 2));

		if (data->interrupt_flag.load()) return;
	}

	{
		uint8_t //top/cur/bot + first/second/last/pre-last
			tf = isCell(-width), ts = isCell(-width + 1),
			tl = isCell(-width + lastElement),
			cf = isCell(0), cs = isCell(1),
			cl = isCell(lastElement);

		for (uint32_t rowIndex_ = 0; rowIndex_ < rowCount; rowIndex_++) {
			const uint32_t rowIndex = rowIndex_ + startRow;
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
		uint8_t //top/cur/bot + first/second/last/pre-last
			tf = isCell(-width),
			tl = isCell(-width + lastElement), tp = isCell(-width + lastElement - 1),
			cf = isCell(0),
			cl = isCell(lastElement), cp = isCell(lastElement - 1);

		for (uint32_t rowIndex_ = 0; rowIndex_ < rowCount; rowIndex_++) {
			const uint32_t rowIndex = rowIndex_ + startRow;
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
	glBufferSubData(GL_SHADER_STORAGE_BUFFER, data->currentOffset() + startIndex, sizeof(FieldCell) * (endIndex - startIndex), &grid->bufferCellAt(startIndex));
	
	//uint8_t a[12];
	//std::fill_n(a, 12, misc::to_underlying<>(FieldCell::WALL));

	//if(data->isOffset)glBufferSubData(GL_SHADER_STORAGE_BUFFER, data->currentOffset() + startIndex, 10, &a);
	
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, bufferP);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

	glFinish();
	data->bufferSend.add(t2.elapsedTime());

	data->grid__iteration++;
	if (data->grid__iteration % (500 + data->index) == 0) std::cout << "grid task " << data->index << ": " << data->gridUpdate.median() << ' ' << data->waiting.median() << ' ' << data->bufferSend.median() << std::endl;

	gpuBufferLock_flag.store(false);

	//bufferSend__.add(t2.elapsedTime());
	//grid__iteration++;
	//if (grid__iteration % 50 == 0) std::cout << grid__.median() << ' ' << bufferSend__.median() << std::endl;
}

//void threadUpdateGrid(std::unique_ptr<Field::GridData>& data) {
//	Timer<> t{};
//	auto& grid = data->grid;
//	const int32_t width = static_cast<int32_t>(grid->width);
//	const int32_t height = static_cast<int32_t>(grid->height);
//	const int lastElement = width - 1;
//	const auto rowCount = data->rowCount;
//	const auto startRow = data->startRow;
//
//	GLFWwindow* window = data->offscreen_context;
//	glfwMakeContextCurrent(window);
//
//	const auto setBufferCellAt = [&grid](uint32_t index, FieldCell cell) -> void { grid->bufferCellAt(index) = cell; }; // { return data.grid[data.gridO.normalizeIndex(index + (data.startRow * data.gridO.width))]; };
//	const auto cellAt = [&grid](int32_t index) -> FieldCell& { return grid->cellAt(index); }; // { return data.grid[data.gridO.normalizeIndex(index + (data.startRow * data.gridO.width))]; };
//	const auto isCell = [&cellAt](int32_t index) -> bool { return fieldCell::isAlive(cellAt(index)); };
//
//	const auto calcNewGenBatch128 = [&grid](uint32_t index) -> __m128i {
//		const auto cellAt = [&grid](uint32_t index) -> FieldCell& { return grid->cellAt(index); }; // { return data.grid[data.gridO.normalizeIndex(index + (data.startRow * data.gridO.width))]; };
//		const auto width = grid->width;
//
//		const __m128i topCellRow = _mm_loadu_si128(reinterpret_cast<__m128i*>(&cellAt(index - width)));
//		const __m128i curCellRow = _mm_loadu_si128(reinterpret_cast<__m128i*>(&cellAt(index)));
//		const __m128i botCellRow = _mm_loadu_si128(reinterpret_cast<__m128i*>(&cellAt(index + width)));
//
//		const __m128i aliveMask = _mm_set1_epi8(0b1);
//		const __m128i lowerForMask = _mm_set1_epi8(0b1111);
//		const __m128i wallCell = _mm_set1_epi8(misc::to_underlying(FieldCell::WALL));
//		const __m128i deadCell = _mm_set1_epi8(misc::to_underlying(FieldCell::DEAD));
//		const __m128i aliveCell = _mm_set1_epi8(misc::to_underlying(FieldCell::ALIVE));
//
//		const __m128i topAliveRow = _mm_and_si128(topCellRow, aliveMask);
//		const __m128i curAliveRow = _mm_and_si128(curCellRow, aliveMask);
//		const __m128i botAliveRow = _mm_and_si128(botCellRow, aliveMask);
//
//		const __m128i cellsAliveInRows = _mm_add_epi8(topAliveRow, _mm_add_epi8(curAliveRow, botAliveRow));
//		const __m128i cells3by3NeighboursAlive =
//			_mm_add_epi8(
//				cellsAliveInRows,
//				_mm_add_epi8(
//					_mm_slli_si128(cellsAliveInRows, 1), // << 8 
//					_mm_srli_si128(cellsAliveInRows, 1)  // >> 8
//				)
//			);
//
//		const __m128i cellsActualNeighboursAlive = _mm_and_si128(_mm_sub_epi8(cells3by3NeighboursAlive, curCellRow), lowerForMask);
//
//		const __m128i walls = _mm_and_si128(_mm_cmpeq_epi8(curAliveRow, wallCell), aliveMask);
//		const __m128i dead = _mm_and_si128(_mm_cmpeq_epi8(curAliveRow, deadCell), aliveMask);
//		const __m128i alive = _mm_and_si128(_mm_cmpeq_epi8(curAliveRow, aliveCell), aliveMask);
//		const __m128i is2Alive = _mm_and_si128(_mm_cmpeq_epi8(cellsActualNeighboursAlive, _mm_set1_epi8(2)), aliveMask);
//		const __m128i is3Alive = _mm_and_si128(_mm_cmpeq_epi8(cellsActualNeighboursAlive, _mm_set1_epi8(3)), aliveMask);
//
//		const __m128i newGen =
//			_mm_or_si128(
//				_mm_slli_epi16(walls, 4),
//				_mm_or_si128(
//					_mm_and_si128(
//						alive,
//						_mm_or_si128(is2Alive, is3Alive)
//					),
//					_mm_and_si128(
//						dead,
//						is3Alive
//					)
//				)
//			);
//
//		return newGen;
//	};
//
//	const uint32_t startIndex = startRow * width;
//	const uint32_t endIndex = (startRow + rowCount) * width;
//
//	const size_t sizeOfBatch = sizeof(__m128i); //dependent on CellsGrid fullSize
//	static_assert(sizeOfBatch > 2, "sizeOfBatch of 0, 1, 2 doesnt make sense here");
//
//	uint32_t i = startIndex;
//
//	for (; (i + sizeOfBatch-2) < endIndex; i+= sizeOfBatch-2) {
//		const __m128i newGen = calcNewGenBatch128(i);
//
//		std::memcpy(&grid->bufferCellAt(i + 1), ((FieldCell*)(&newGen)) + 1, sizeof(FieldCell) * (sizeOfBatch - 2));
//
//		if (data->interrupt_flag.load()) return;
//	}
//
//	if (i != endIndex) {
//		uint32_t remainingCells = endIndex - i;
//
//		const __m128i newGen = calcNewGenBatch128(i);
//
//		std::memcpy(&grid->bufferCellAt(i + 1), ((FieldCell*)(&newGen)) + 1, sizeof(FieldCell) * (remainingCells - 2));
//
//		if (data->interrupt_flag.load()) return;
//	}
//
//	{
//		uint8_t //top/cur/bot + first/second/last/pre-last
//			tf = isCell(-width), ts = isCell(-width + 1),
//			tl = isCell(-width + lastElement),
//			cf = isCell(0), cs = isCell(1),
//			cl = isCell(lastElement);
//
//		for (uint32_t rowIndex_ = 0; rowIndex_ < rowCount; rowIndex_++) {
//			const uint32_t rowIndex = rowIndex_ + startRow;
//			const uint32_t row = rowIndex * width;
//
//			uint8_t
//				bf = isCell(row + width),
//				bs = isCell(row + width + 1),
//				bl = isCell(row + width + lastElement);
//			//first element
//			{
//				const int index = row;
//				const auto curCell = cellAt(index);
//				const unsigned char    topRowNeighbours = tf + ts + tl;
//				const unsigned char bottomRowNeighbours = bf + bs + bl;
//				const unsigned char aliveNeighbours = topRowNeighbours + cl + cs + bottomRowNeighbours;
//
//				setBufferCellAt(index, fieldCell::nextGeneration(curCell, aliveNeighbours));
//			}
//
//			tf = cf;
//			ts = cs;
//			tl = cl;
//			cf = bf;
//			cs = bs;
//			cl = bl;
//		}
//	}
//	if (data->interrupt_flag.load()) return;
//
//	{
//		uint8_t //top/cur/bot + first/second/last/pre-last
//			tf = isCell(-width),
//			tl = isCell(-width + lastElement), tp = isCell(-width + lastElement - 1),
//			cf = isCell(0),
//			cl = isCell(lastElement), cp = isCell(lastElement - 1);
//
//		for (uint32_t rowIndex_ = 0; rowIndex_ < rowCount; rowIndex_++) {
//			const uint32_t rowIndex = rowIndex_ + startRow;
//			const uint32_t row = rowIndex * width;
//
//			uint8_t
//				bf = isCell(row + width),
//				bp = isCell(row + width + lastElement - 1),
//				bl = isCell(row + width + lastElement);
//
//			{ //last element
//				const int index = row + lastElement;
//				const auto curCell = cellAt(index);
//				const unsigned char    topRowNeighbours = tp + tl + tf;
//				const unsigned char bottomRowNeighbours = bp + bl + bf;
//				const unsigned char aliveNeighbours = topRowNeighbours + cp + cf + bottomRowNeighbours;
//
//				setBufferCellAt(index, fieldCell::nextGeneration(curCell, aliveNeighbours));
//			}
//
//			tf = cf;
//			tp = cp;
//			tl = cl;
//			cf = bf;
//			cp = bp;
//			cl = bl;
//		}
//	}
//	if (data->interrupt_flag.load()) return;
//
//	data->gridUpdate.add(t.elapsedTime());
//
//
//	Timer<> t2{};
//	auto &gpuBufferLock_flag = data->gpuBufferLock_flag;
//	bool expected = false;
//	while (!gpuBufferLock_flag.compare_exchange_weak(expected, true)) {}// ::std::cout << data->index << ":waiting" << ::std::endl; }
//
//	data->waiting.add(t2.elapsedTime());
//	
//	const auto bufferP = data->bufferP;
//	glBindBuffer(GL_SHADER_STORAGE_BUFFER, bufferP);
//	glBufferSubData(GL_SHADER_STORAGE_BUFFER, data->currentOffset() + startIndex, sizeof(FieldCell) * (endIndex - startIndex), &grid->bufferCellAt(startIndex));
//	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, bufferP);
//	glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
//
//	glFinish();
//	data->bufferSend.add(t2.elapsedTime());
//
//	data->grid__iteration++;
//	if (data->grid__iteration % (200 + data->index) == 0) std::cout << data->index << ":" << data->gridUpdate.median() << ' ' << data->waiting.median() << ' ' << data->bufferSend.median() << std::endl;
//
//	gpuBufferLock_flag.store(false);
//
//	//bufferSend__.add(t2.elapsedTime());
//	//grid__iteration++;
//	//if (grid__iteration % 50 == 0) std::cout << grid__.median() << ' ' << bufferSend__.median() << std::endl;
//}

FieldCell* Field::grid() {
	shouldUpdateGrid = false;
	return &this->gridPimpl->cellAt(0);
}

Field::Field(const uint32_t gridWidth, const uint32_t gridHeight, const size_t numberOfTasks_, const GLuint bufferP_, GLFWwindow* window) :
	gridPimpl{ new FieldPimpl(gridWidth, gridHeight) },
	shouldUpdateGrid(true),
	numberOfTasks(numberOfTasks_),
	gridTasks{ new std::unique_ptr<Task<std::unique_ptr<GridData>>>[numberOfTasks_] },
	indecesToBrokenCells{ },
	bufferP{ bufferP_ }
{
	assert(numberOfTasks >= 1);

	//grid__.median() == 35 ms, bufferSend__.median() == 110 ms.In my case
	//const double fraction = 35.0 / 110.0; // fractionOfWorkTimeToSendData

	const double fraction = 1; // fractionOfWorkTimeToSendData
	/*
	find amountOfWorkT1. given fraction, workload (grid_height)
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
	*/
	const double amountOfWorkT1 = height() / ((1.0 + fraction) * (pow((1.0 + fraction), numberOfTasks - 1) - 1)) * fraction;

	const auto createGridTask = [window, this](const uint32_t index, const uint32_t rowsBefore, const uint32_t numberOfRows) -> void {
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
						rowsBefore,
						numberOfRows,
						bufferP,
						isFieldGPUBufferOffset,
						size(),
						gpuBufferLock_flag,
						taskWindow
					}
				)
			}
		);
	};

	uint32_t rowsBefore = 0;
	uint32_t remainingRows = height();

	/*for (uint32_t i = 0; i < numberOfTasks; i++) {
		const uint32_t numberOfRows = remainingRows / (numberOfTasks - i);

		createGridTask(i, rowsBefore, numberOfRows);

		rowsBefore += numberOfRows;
		remainingRows -= numberOfRows;
	}*/

	double currentTaskWork = amountOfWorkT1;

	for (uint32_t i = 0; i < numberOfTasks - 1; i++) {
		const uint32_t numberOfRows = misc::min(remainingRows, static_cast<uint32_t>(ceil(currentTaskWork)));
		::std::cout << numberOfRows << std::endl;

		createGridTask(i, rowsBefore, numberOfRows);

		rowsBefore += numberOfRows;
		remainingRows -= numberOfRows;

		currentTaskWork *= (1.0 + fraction);
	}

	{//last
		const uint32_t numberOfRows = remainingRows;
		::std::cout << numberOfRows << std::endl;

		createGridTask(numberOfTasks - 1, rowsBefore, numberOfRows);

		rowsBefore += numberOfRows;
		remainingRows -= numberOfRows;
	}

	assert(remainingRows == 0);
	assert(rowsBefore == height());

	deployGridTasks();

}
Field::~Field() = default;

uint32_t Field::currentOffset() {
	return size() * isFieldGPUBufferOffset;
}

uint32_t Field::nextOffset() {
	return size() * !isFieldGPUBufferOffset;
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
	glBufferSubData(GL_SHADER_STORAGE_BUFFER, nextOffset() + 0, sizeof(FieldCell) * size(), &gridPimpl->cellAt(0));
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
		glBufferSubData(GL_SHADER_STORAGE_BUFFER, nextOffset() + index, sizeof(FieldCell), &gridPimpl->cellAt(index));
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
					glBufferSubData(GL_SHADER_STORAGE_BUFFER, currentOffset() + offsetedIndex, sizeof(FieldCell), &gridPimpl->bufferCellAt(offsetedIndex));
				}
			}

		}
		indecesToBrokenCells.clear();

		glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, bufferP);
		glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
	}

	gridPimpl->prepareNext();
	shouldUpdateGrid = true;
	isFieldGPUBufferOffset = !isFieldGPUBufferOffset;
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