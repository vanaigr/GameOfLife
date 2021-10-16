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
		startRow(startRow_),
		rowCount(rowCount_),
		 bufferP(bufferP_),
		 isOffset(isOffset_),
		 offset(offset_),
		gpuBufferLock_flag(gpuBufferLock_flag_),
		offscreen_context(offscreen_context_)
	{}
};

//struct Field::PackedGridData {
//	std::unique_ptr<FieldPimpl>& grid;
//	std::shared_ptr<uint32_t[]>& packedGrid;
//
//	uint32_t startCell;
//	uint32_t cellsCount;
//};

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

void threadUpdateGrid(std::unique_ptr<Field::GridData>& data) {
	Timer<> t{};
	auto& grid = data->grid;
	const int32_t width = grid->width;
	const int32_t height = grid->height;
	const int lastElement = width - 1;
	const uint32_t rowCount = data->rowCount;
	const uint32_t startRow = data->startRow;

	GLFWwindow* window = data->offscreen_context;
	glfwMakeContextCurrent(window);

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

	data->gridUpdate.add(t.elapsedTime());


	Timer<> t2{};
	auto &gpuBufferLock_flag = data->gpuBufferLock_flag;
	bool expected = false;
	while (!gpuBufferLock_flag.compare_exchange_weak(expected, true)) {}// ::std::cout << data->index << ":waiting" << ::std::endl; }

	data->waiting.add(t2.elapsedTime());
	
	const auto bufferP = data->bufferP;
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, bufferP);
	glBufferSubData(GL_SHADER_STORAGE_BUFFER, data->currentOffset() + startIndex, sizeof(FieldCell) * (endIndex - startIndex), &grid->bufferCellAt(startIndex));
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, bufferP);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

	glFinish();
	data->bufferSend.add(t2.elapsedTime());

	data->grid__iteration++;
	if (data->grid__iteration % (200 + data->index) == 0) std::cout << data->index << ":" << data->gridUpdate.median() << ' ' << data->waiting.median() << ' ' << data->bufferSend.median() << std::endl;

	gpuBufferLock_flag.store(false);

	//bufferSend__.add(t2.elapsedTime());
	//grid__iteration++;
	//if (grid__iteration % 50 == 0) std::cout << grid__.median() << ' ' << bufferSend__.median() << std::endl;
}

FieldCell* Field::grid() {
	shouldUpdateGrid = false;
	return &this->gridPimpl->cellAt(0);
}

Field::Field(const uint32_t gridWidth, const uint32_t gridHeight, const size_t numberOfTasks_, const GLuint bufferP, GLFWwindow* window) :
	gridPimpl{ new FieldPimpl(gridWidth, gridHeight) },
	//packedGridSize(size() / (32 / 2) + 1),
	shouldUpdateGrid(false),

	numberOfTasks(numberOfTasks_),
	gridTasks{ new std::unique_ptr<Task<std::unique_ptr<GridData>>>[numberOfTasks_] },
	indecesToBrokenCells{ }
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

	const auto createGridTask = [window, this, &bufferP](const uint32_t index, const uint32_t rowsBefore, const uint32_t numberOfRows) -> void {
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