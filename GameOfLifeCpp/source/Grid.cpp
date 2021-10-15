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



static const bool isCellAlive(uint32_t index, const Field &grid)  {
	return grid.cellAtIndex(index) == FieldCell::ALIVE; }
;

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
		fullSize{ size + padding + padding + 16u },//+2 duplicate rows (top row duplicate before first row, bottom row duplicate after last row) + extra bytes for memcopy/vector instruction in thread grid update
		current{ new FieldCell[fullSize * 2]{ FieldCell::DEAD } }, 
		buffer{ new FieldCell[fullSize * 2 /*tmp '*2' */]{ FieldCell::DEAD } }
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
	std::unique_ptr<FieldPimpl>& grid;
	std::atomic_bool& interrupt_flag;
	uint32_t startRow;
	uint32_t rowCount;

public:
	GridData(
		std::unique_ptr<FieldPimpl>& grid_,
		std::atomic_bool& interrupt_flag_,
		uint32_t startRow_,
		uint32_t rowCount_
	) :
		grid(grid_),
		interrupt_flag(interrupt_flag_),
		startRow(startRow_),
		rowCount(rowCount_)
	{}
};

struct Field::PackedGridData {
	std::unique_ptr<FieldPimpl>& grid;
	std::shared_ptr<uint32_t[]>& packedGrid;

	uint32_t startCell;
	uint32_t cellsCount;
};

FieldCell updatedCell_(const int32_t index, const std::unique_ptr<Field::FieldPimpl>& cellsGrid) {
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

		if (cell == FieldCell::DEAD && aliveNeighbours == 3) {
			return FieldCell::ALIVE;
		}
		else if (cell == FieldCell::ALIVE && (aliveNeighbours < 2 || aliveNeighbours > 3)) {
			return FieldCell::DEAD;
		}
	}

	return cell;
}

MedianCounter grid__{};
int grid__iteration = 0;

void threadUpdateGrid(std::unique_ptr<Field::GridData>& data) {
	
	Timer<> t{};
	auto& grid = data->grid;
	const int32_t width = grid->width;
	const int32_t height = grid->height;
	const int lastElement = width - 1;
	const int rowCount = data->rowCount;
	const int startRow = data->startRow;

	const auto setBufferCellAt = [&grid](uint32_t index, FieldCell cell) -> void { grid->bufferCellAt(index) = cell; }; // { return data.grid[data.gridO.normalizeIndex(index + (data.startRow * data.gridO.width))]; };
	const auto cellAt = [&grid](int32_t index) -> FieldCell& { return grid->cellAt(index); }; // { return data.grid[data.gridO.normalizeIndex(index + (data.startRow * data.gridO.width))]; };
	const auto isCell = [&cellAt](int32_t index) -> bool { return fieldCell::isAlive(cellAt(index)); };

	const auto calcNewGenBatch128 = [&grid](uint32_t index) -> __m128i {
		const auto cellAt = [&grid](uint32_t index) -> FieldCell& { return grid->cellAt(index); }; // { return data.grid[data.gridO.normalizeIndex(index + (data.startRow * data.gridO.width))]; };
		const auto width = grid->width;

		const __m128i topCellRow = _mm_loadu_si128(reinterpret_cast<__m128i*>(&cellAt(index - width)));
		const __m128i curCellRow = _mm_loadu_si128(reinterpret_cast<__m128i*>(&cellAt(index)));
		const __m128i botCellRow = _mm_loadu_si128(reinterpret_cast<__m128i*>(&cellAt(index + width)));

		const __m128i aliveMask = _mm_set1_epi8(0b1);
		const __m128i lowerForMask = _mm_set1_epi8(0b1111);
		const __m128i two = _mm_set1_epi8(2);
		const __m128i three = _mm_set1_epi8(3);
		const __m128i wallCell = _mm_set1_epi8(misc::to_underlying(FieldCell::WALL));
		const __m128i deadCell = _mm_set1_epi8(misc::to_underlying(FieldCell::DEAD));
		const __m128i aliveCell = _mm_set1_epi8(misc::to_underlying(FieldCell::ALIVE));

		const __m128i topAliveRow = _mm_and_si128(topCellRow, aliveMask);
		const __m128i curAliveRow = _mm_and_si128(curCellRow, aliveMask);
		const __m128i botAliveRow = _mm_and_si128(botCellRow, aliveMask);

		const __m128i cellsAliveInRows = _mm_add_epi8(topAliveRow, _mm_add_epi8(curAliveRow, botAliveRow));
		const __m128i cells3by3NeighboursAlive =
			_mm_add_epi8(
				cellsAliveInRows,
				_mm_add_epi8(
					_mm_slli_si128(cellsAliveInRows, 1), // << 8 
					_mm_srli_si128(cellsAliveInRows, 1)  // >> 8
				)
			);

		const __m128i cellsActualNeighboursAlive = _mm_and_si128(_mm_sub_epi8(cells3by3NeighboursAlive, curCellRow), lowerForMask);

		const __m128i walls = _mm_and_si128(_mm_cmpeq_epi8(curAliveRow, wallCell), aliveMask);
		const __m128i dead = _mm_and_si128(_mm_cmpeq_epi8(curAliveRow, deadCell), aliveMask);
		const __m128i alive = _mm_and_si128(_mm_cmpeq_epi8(curAliveRow, aliveCell), aliveMask);
		const __m128i is2Alive = _mm_and_si128(_mm_cmpeq_epi8(cellsActualNeighboursAlive, two), aliveMask);
		const __m128i is3Alive = _mm_and_si128(_mm_cmpeq_epi8(cellsActualNeighboursAlive, three), aliveMask);

		const __m128i newGen =
			_mm_or_si128(
				_mm_slli_epi16(walls, 4),
				_mm_or_si128(
					_mm_and_si128(
						alive,
						_mm_or_si128(is2Alive, is3Alive)
					),
					_mm_and_si128(
						dead,
						is3Alive
					)
				)
			);

		return newGen;
	};

	const uint32_t startIndex = startRow * width;
	const uint32_t endIndex = (startRow + rowCount) * width;

	const size_t sizeOfBatch = sizeof(__m128i); //dependent on CellsGrid fullSize
	static_assert(sizeOfBatch > 2, "sizeOfBatch of 0, 1, 2 doesnt make sense here");

	uint32_t i = startIndex;

	for (; (i + sizeOfBatch-2) < endIndex; i+= sizeOfBatch-2) {
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

	/*
	{
		uint8_t //top/cur/bot + first/second/last/pre-last
			tf = isCell(-width), ts = isCell(-width + 1),
			tl = isCell(-width + lastElement), tp = isCell(-width + lastElement - 1),
			cf = isCell(0), cs = isCell(1),
			cl = isCell(lastElement), cp = isCell(lastElement - 1);

		for (uint32_t rowIndex_ = 0; rowIndex_ < rowCount; rowIndex_++) {
			const uint32_t rowIndex = rowIndex_ + startRow;
			const uint32_t row = rowIndex * width;

			uint8_t
				bf = isCell(row + width),
				bs = isCell(row + width + 1),
				bp = isCell(row + width + lastElement - 1),
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
			{ //last element
				const int index = row + lastElement;
				const auto curCell = cellAt(index);
				const unsigned char    topRowNeighbours = tp + tl + tf;
				const unsigned char bottomRowNeighbours = bp + bl + bf;
				const unsigned char aliveNeighbours = topRowNeighbours + cp + cf + bottomRowNeighbours;

				setBufferCellAt(index, fieldCell::nextGeneration(curCell, aliveNeighbours));
			}
			if (data->interrupt_flag.load()) return;

			tf = cf;
			ts = cs;
			tp = cp;
			tl = cl;
			cf = bf;
			cs = bs;
			cp = bp;
			cl = bl;
		}
	}
	/*
	for (uint32_t rowIndex_ = 0; rowIndex_ < rowCount; rowIndex_++) {
		const uint32_t rowIndex = rowIndex_ + startRow;
		const uint32_t row = rowIndex * width;

		//first element
		{
			const int index = row;
			const auto curCell = cellAt(index);
			const unsigned char bottomRowNeighbours = isCell(row + width + lastElement) + isCell(index + width) + isCell(index + width + 1);
			const unsigned char    topRowNeighbours = isCell(row - width + lastElement) + isCell(index - width) + isCell(index - width + 1);
			const unsigned char aliveNeighbours = topRowNeighbours + isCell(row + lastElement) + isCell(index + 1) + bottomRowNeighbours;

			setBufferCellAt(index, fieldCell::nextGeneration(curCell, aliveNeighbours));
		}
		{ //last element
			const int index = row + lastElement;
			const auto curCell = cellAt(index);
			const unsigned char bottomRowNeighbours = isCell(index + width - 1) + isCell(index + width) + isCell(row + width + 0);
			const unsigned char    topRowNeighbours = isCell(index - width - 1) + isCell(index - width) + isCell(row - width + 0);
			const unsigned char aliveNeighbours = topRowNeighbours + isCell(index - 1) + isCell(row + 0) + bottomRowNeighbours;

			setBufferCellAt(index, fieldCell::nextGeneration(curCell, aliveNeighbours));
		}
		if (data->interrupt_flag.load()) return;
	}*/

	//glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssboHandle());
	//glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, endIndex - startIndex, &cellAt(startIndex));
	//glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, ssboHandle());
	//glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

	grid__.add(t.elapsedTime());
	grid__iteration++;
	//if (grid__iteration % 50 == 0) std::cout << grid__.median() << std::endl;
}

/*void threadUpdateGrid(std::unique_ptr<Field::GridData>& data) {
	Timer<> t{};
	auto& grid = data->grid;
	const auto width = grid->width;
	const auto height = grid->height;
	const int lastElement = width - 1;
	const int rowCount = data->rowCount;
	const int startRow = data->startRow;

	const auto setBufferCellAt = [&grid](uint32_t index, FieldCell cell) -> void { grid->bufferCellAt(index) = cell; }; // { return data.grid[data.gridO.normalizeIndex(index + (data.startRow * data.gridO.width))]; };
	const auto cellAt = [&grid](uint32_t index) -> FieldCell& { return grid->cellAt(index); }; // { return data.grid[data.gridO.normalizeIndex(index + (data.startRow * data.gridO.width))]; };
	const auto isCell = [&cellAt](uint32_t index) -> bool { return fieldCell::isAlive(cellAt(index)); };

	for (uint32_t rowIndex_ = 0; rowIndex_ < rowCount; rowIndex_++) {
		const uint32_t rowIndex = rowIndex_ + startRow;
		const uint32_t row = rowIndex * width;

		uint32_t startColumn = 0;
		for (;;) {
			const uint32_t startIndex = row + startColumn;

			const size_t sizeOfBatch = 64; //dependent on CellsGrid fullSize
			static_assert(sizeOfBatch > 2, "sizeOfBatch of 0, 1, 2 doesnt make sense here");
			FieldCell topCellRow[sizeOfBatch];
			FieldCell curCellRow[sizeOfBatch];
			FieldCell botCellRow[sizeOfBatch];

			std::memcpy(&topCellRow, &cellAt(startIndex - width), sizeof(topCellRow));
			std::memcpy(&curCellRow, &cellAt(startIndex), sizeof(curCellRow));
			std::memcpy(&botCellRow, &cellAt(startIndex + width), sizeof(botCellRow));

			uint8_t cellsAliveInRows[sizeOfBatch];
			for (size_t i = 0; i < sizeOfBatch; i++) cellsAliveInRows[i] = fieldCell::isAlive(topCellRow[i]) + fieldCell::isAlive(curCellRow[i]) + fieldCell::isAlive(botCellRow[i]);

			uint8_t cellsNeighboursAlive[sizeOfBatch-2];
			for (size_t i = 1; i < sizeOfBatch-1; i++) cellsNeighboursAlive[i-1] = cellsAliveInRows[i-1] + cellsAliveInRows[i] + cellsAliveInRows[i + 1];

			uint32_t remainingCells = width - startColumn;
			uint32_t maxCells = misc::min(remainingCells, sizeOfBatch);

			FieldCell newGen[sizeOfBatch-2];
			for (uint32_t offset = 0; offset < maxCells - 2; offset++) {
				const uint32_t cellOffset = offset + 1;
				const uint32_t curCellIndex = startIndex + cellOffset;
				const auto curCell = cellAt(curCellIndex);
				const uint8_t aliveNeighbours = cellsNeighboursAlive[offset] - fieldCell::isAlive(curCellRow[cellOffset]);

				//setBufferCellAt(curCellIndex, gridCell::nextGeneration(curCell, aliveNeighbours));
				newGen[offset] = fieldCell::nextGeneration(curCell, aliveNeighbours);
			}

			std::memcpy(&grid->bufferCellAt(startIndex + 1), &newGen, sizeof(FieldCell) * (maxCells - 2));

			if (remainingCells < sizeOfBatch) {
				break;
			}

			startColumn += sizeOfBatch-2;
		}

		//first element
		{
			const int index = row;
			const auto curCell = cellAt(index);
			const unsigned char bottomRowNeighbours = isCell(row + width + lastElement) + isCell(index + width) + isCell(index + width + 1);
			const unsigned char    topRowNeighbours = isCell(row - width + lastElement) + isCell(index - width) + isCell(index - width + 1);
			const unsigned char aliveNeighbours = topRowNeighbours + isCell(row + lastElement) + isCell(index + 1) + bottomRowNeighbours;

			setBufferCellAt(index, fieldCell::nextGeneration(curCell, aliveNeighbours));
		}
		{ //last element
			const int index = row + lastElement;
			const auto curCell = cellAt(index);
			const unsigned char bottomRowNeighbours = isCell(index + width - 1) + isCell(index + width) + isCell(row + width + 0);
			const unsigned char    topRowNeighbours = isCell(index - width - 1) + isCell(index - width) + isCell(row - width + 0);
			const unsigned char aliveNeighbours = topRowNeighbours + isCell(index - 1) + isCell(row + 0) + bottomRowNeighbours;

			setBufferCellAt(index, fieldCell::nextGeneration(curCell, aliveNeighbours));
		}
		if (data->interrupt_flag.load()) return;
	}

	grid__.add(t.elapsedTime());
	grid__iteration++;
	if (grid__iteration % 50 == 0) std::cout << grid__.median() << std::endl;
}*/

FieldCell* Field::grid() {
	shouldUpdateGrid = false;
	return &this->gridPimpl->cellAt(0);
}

Field::Field(const uint32_t gridWidth, const uint32_t gridHeight, const size_t numberOfTasks_) :
	gridPimpl{ new FieldPimpl(gridWidth, gridHeight) },
	packedGridSize((gridWidth * gridHeight) / (32 / 2) + 1),
	shouldUpdateGrid(false),

	numberOfTasks(numberOfTasks_),
	gridTasks{ NULL },
	indecesToBrokenCells{ }
{
	assert(numberOfTasks >= 1);

	gridTasks.reset(new std::unique_ptr<Task<std::unique_ptr<GridData>>>[numberOfTasks]);
	uint32_t rowsBefore = 0;
	uint32_t remainingRows = height();
	for (uint32_t i = 0; i < numberOfTasks; i++) {
		const uint32_t numberOfRows = remainingRows / (numberOfTasks - i);

		gridTasks.get()[i] = std::unique_ptr<Task<std::unique_ptr<GridData>>>(
			new Task<std::unique_ptr<GridData>>{
				threadUpdateGrid,
				std::make_unique<GridData>(
					GridData{
						this->gridPimpl,
						this->interrupt_flag,
						rowsBefore,
						numberOfRows
					}
				)
			}

		);

		rowsBefore += numberOfRows;
		remainingRows -= numberOfRows;
	}

	deployGridTasks();

}
Field::~Field() = default;


void Field::fill(const FieldCell cell) {
	interrupt_flag.store(true);
	waitForGridTasks();
	interrupt_flag.store(false);

	gridPimpl->fill(cell);

	shouldUpdateGrid = true;

	indecesToBrokenCells.clear();
	deployGridTasks();
}

FieldCell Field::cellAtIndex(const uint32_t index) const {
	return gridPimpl->cellAt(index);
}

void Field::setCellAtIndex(const uint32_t index, FieldCell cell) {
	auto& gridCell = gridPimpl->cellAt(index);
	if (gridCell != cell) {
		gridCell = cell;
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
		for (auto it = indecesToBrokenCells.begin(); it != indecesToBrokenCells.end(); it++) {
			const uint32_t index = *it;

			for (int yo = -1; yo <= 1; yo++) {
				for (int xo = -1; xo <= 1; xo++) {
					int x = ((index + xo) + width()) % width(),
						y = (((index / width()) + yo) + height()) % height();
					int offsetedIndex = x + width() * y;
					const auto newCell = updatedCell_(offsetedIndex, gridPimpl);
					gridPimpl->bufferCellAt(offsetedIndex) = newCell;
				}
			}

		}
		indecesToBrokenCells.clear();
	}

	gridPimpl->prepareNext();
	shouldUpdateGrid = true;

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