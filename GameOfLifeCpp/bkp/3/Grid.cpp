#include "Misc.h"
#include "Grid.h"
#include <vector>

#include <cassert> 


static const auto updatedCell = [](GridCell cell, unsigned char aliveNeighbours) {
	if (cell == GridCell::DEAD && aliveNeighbours == 3) {
		return GridCell::ALIVE;
	}
	else if (cell == GridCell::ALIVE && (aliveNeighbours < 2 || aliveNeighbours > 3)) {
		return GridCell::DEAD;
	}
	return cell;
};

static const bool isCellAlive(uint32_t index, const Grid &grid)  {
	return grid.cellAtIndex(index) == GridCell::ALIVE; }
;

template<uint32_t rowStartOffset, uint32_t rowEndOffset>
static inline void updateGrid_noVerticalBoundaries(TaskData& data) {
	const auto& gridO = data.gridO;
	//const auto& grid = data.grid;
	const &gl = data.gridLocal;
	const auto width = gridO.width;
	const auto height = gridO.height;
	const int lastElement = width - 1;
	const int lastRow = width * (height - 1);
	const auto isCell = [&gridO](uint32_t index) -> bool { return gl };

	for (uint32_t rowOffset = 0 + rowStartOffset; rowOffset < data.rowCount - rowEndOffset; rowOffset++) {
		const uint32_t outputRow = rowOffset * width;
		const uint32_t row = data.startRow + rowOffset;
		const uint32_t rowIndex = row * width;

		for (int col = 1; col < width - 1; col++) {
			const uint32_t outputIndex = outputRow + col;
			const int index = col + rowIndex;
			const auto curCell = grid[index];
			if (curCell != GridCell::WALL) {
				const unsigned char bottomRowNeighbours = isCell(index + width - 1) + isCell(index + width) + isCell(index + width + 1);
				const unsigned char topRowNeighbours = isCell(index - width - 1) + isCell(index - width) + isCell(index - width + 1);
				const unsigned char aliveNeighbours = topRowNeighbours + isCell(index - 1) + isCell(index + 1) + bottomRowNeighbours;

				data.output[outputIndex] = updatedCell(curCell, aliveNeighbours);
			}
		}

		//first element
		{
			const int index = rowIndex;
			const auto curCell = grid[index];
			if (curCell != GridCell::WALL) {
				const unsigned char bottomRowNeighbours = isCell(index + width + lastElement) + isCell(index + width) + isCell(index + width + 1);
				const unsigned char topRowNeighbours = isCell(index - width + lastElement) + isCell(index - width) + isCell(index - width + 1);
				const unsigned char aliveNeighbours = topRowNeighbours + isCell(index + lastElement) + isCell(index + 1) + bottomRowNeighbours;

				data.output[outputRow] = updatedCell(curCell, aliveNeighbours);
			}
		}
		{ //last element
			const int index = rowIndex + lastElement;
			const auto curCell = grid[index];
			if (curCell != GridCell::WALL) {
				const unsigned char bottomRowNeighbours = isCell(index + width - 1) + isCell(index + width) + isCell(rowIndex + width);
				const unsigned char    topRowNeighbours = isCell(index - width - 1) + isCell(index - width) + isCell(rowIndex - width);
				const unsigned char aliveNeighbours = topRowNeighbours + isCell(index - 1) + isCell(rowIndex + 0) + bottomRowNeighbours;

				data.output[outputRow + lastElement] = updatedCell(curCell, aliveNeighbours);

			}
		}
	}
}

void threadUpdateGridTop(TaskData& data) {
	updateGrid_noVerticalBoundaries<1, 0>(data);

	const auto& gridO = data.gridO;
	const auto& grid = data.grid;
	const auto width = gridO.width;
	const auto height = gridO.height;
	const int lastElement = width - 1;
	const int lastRow = width * (height - 1);
	const auto isCell = [&gridO](uint32_t index) -> bool { return isCellAlive(index, gridO); };

	{//first row
		for (int col = 1; col < width - 1; col++) {
			int index = col;
			const auto curCell = grid[index];
			if (curCell != GridCell::WALL) {
				const unsigned char topRowNeighbours = isCell(lastRow + index - 1) + isCell(lastRow + index) + isCell(lastRow + index + 1);
				const unsigned char bottomRowNeighbours = isCell(width + index - 1) + isCell(index + width) + isCell(width + index + 1);
				const unsigned char aliveNeighbours = topRowNeighbours + isCell(index - 1) + isCell(index + 1) + bottomRowNeighbours;

				data.output[index] = updatedCell(curCell, aliveNeighbours);
			}
		}

		{ //first element
			const auto curCell = grid[0];
			if (curCell != GridCell::WALL) {
				const unsigned char    topRowNeighbours = isCell(lastRow + lastElement) + isCell(lastRow + 0) + isCell(lastRow + 1);
				const unsigned char bottomRowNeighbours = isCell(width + lastElement) + isCell(width) + isCell(width + 1);
				const unsigned char aliveNeighbours = topRowNeighbours + isCell(lastElement) + isCell(1) + bottomRowNeighbours;

				data.output[0] = updatedCell(curCell, aliveNeighbours);
			}
		}
		{ //last element
			int index = lastElement;
			const auto curCell = grid[index];
			if (curCell != GridCell::WALL) {
				const unsigned char bottomRowNeighbours = isCell(lastRow + lastElement - 1) + isCell(lastRow + lastElement) + isCell(lastRow);
				const unsigned char    topRowNeighbours = isCell(width + lastElement - 1) + isCell(width + lastElement) + isCell(width);
				const unsigned char aliveNeighbours = topRowNeighbours + isCell(lastElement - 1) + isCell(0) + bottomRowNeighbours;

				data.output[index] = updatedCell(curCell, aliveNeighbours);
			}
		}
	}
}

void threadUpdateGrid(TaskData& data) {
	updateGrid_noVerticalBoundaries<0, 0>(data);

}
void threadUpdateGridBot(TaskData& data) {
	updateGrid_noVerticalBoundaries<0, 1>(data);


	const auto& gridO = data.gridO;
	const auto& grid = data.grid;
	const auto width = gridO.width;
	const auto height = gridO.height;
	const int lastElement = width - 1;
	const int lastRow = width * (height - 1);
	const auto isCell = [&gridO](uint32_t index) -> bool { return isCellAlive(index, gridO); };

	{//last row
		const uint32_t outputRow = (data.rowCount - 1) * width;
		for (int col = 1; col < width - 1; col++) {
			int index = lastRow + col;
			const uint32_t outputIndex = outputRow + col;
			const auto curCell = grid[index];
			if (curCell != GridCell::WALL) {
				const unsigned char topRowNeighbours = isCell(index - width - 1) + isCell(index - width) + isCell(index - width + 1);
				const unsigned char bottomRowNeighbours = isCell(col - 1) + isCell(col) + isCell(col + 1);
				const unsigned char aliveNeighbours = topRowNeighbours + isCell(lastRow + col - 1) + isCell(lastRow + col + 1) + bottomRowNeighbours;

				data.output[outputIndex] = updatedCell(curCell, aliveNeighbours);
			}
		}

		{ //first element
			int index = lastRow;
			const auto curCell = grid[index];
			if (curCell != GridCell::WALL) {
				const unsigned char    topRowNeighbours = isCell(index - width + lastElement) + isCell(index - width) + isCell(index - width + 1);
				const unsigned char bottomRowNeighbours = isCell(lastElement) + isCell(0) + isCell(1);
				const unsigned char aliveNeighbours = topRowNeighbours + isCell(index + lastElement) + isCell(index + 1) + bottomRowNeighbours;

				data.output[outputRow] = updatedCell(curCell, aliveNeighbours);
			}
		}
		{ //last element
			int index = lastRow + lastElement;
			const auto curCell = grid[index];
			if (curCell != GridCell::WALL) {
				const unsigned char    topRowNeighbours = isCell(index - width - 1) + isCell(index - width) + isCell(lastRow - width);
				const unsigned char bottomRowNeighbours = isCell(lastElement - 1) + isCell(lastElement) + isCell(0);
				const unsigned char aliveNeighbours = topRowNeighbours + isCell(index - 1) + isCell(lastRow) + bottomRowNeighbours;

				data.output[outputRow + lastElement] = updatedCell(curCell, aliveNeighbours);
			}
		}
	}
}

Grid::Grid(const uint32_t gridWidth, const uint32_t gridHeight, const size_t numberOfTasks_) :
	width(gridWidth), height(gridHeight), 

	size(gridWidth* gridHeight), 
	packedGridSize((width* height) / (32 / 2) + 1),
	shouldUpdatePackedGrid(false),

	grid{ new GridCell[size]{ GridCell::DEAD } },
	packedGrid_{ new uint32_t[size]{ 0 } },

	numberOfTasks(numberOfTasks_),
	tasks{ NULL }
{
	assert(numberOfTasks >= 2); //0, 1 thread mod is not supproted

	/*const uint32_t threadGridHeight = gridHeight / numberOfTasks;
	for (uint32_t i = 0; i < numberOfTasks - 1; i++) {
		const uint32_t length = gridWidth * threadGridHeight;
		const uint32_t startIndex = length * i;
		//const auto &output = (threadsData.get()[i] = std::move<>(std::unique_ptr<GridCell[]>{new GridCell[length]{ GridCell::DEAD } }));
		
		tasks[i] = std::thread(threadUpdateGrid, std::ref(*this), std::ref(output), startIndex, length);
	}
	const uint32_t totalSize = (numberOfThreads - 1) * threadGridHeight;
	const uint32_t remainingSize = height - totalSize;
	const uint32_t length = gridWidth * remainingSize;
	const uint32_t startIndex = size - length;
	const auto& output = (threadsData.get()[numberOfThreads - 1] = std::move<>(std::unique_ptr<GridCell[]>{new GridCell[length]{ GridCell::DEAD } }));
	threads[numberOfThreads - 1] = std::thread(threadUpdateGrid, std::ref(*this), std::ref(output), startIndex, length);*/
}

const std::shared_ptr<const uint32_t[]> Grid::packedGrid() {
	if (shouldUpdatePackedGrid) {
		for (size_t i = 0; i < size; i++) {
			const auto cellAtIndex_ = grid.get()[i];
			unsigned int isWall = cellAtIndex_ == GridCell::WALL;
			unsigned int isCell = cellAtIndex_ == GridCell::ALIVE;
			unsigned int mask = (((isWall << 1) & 0b10) | isCell) & 0b11;

			unsigned int arrIndex = i / (32 / 2);
			unsigned int arrShift = (i % (32 / 2)) * 2;
			packedGrid_[arrIndex] = (packedGrid_[arrIndex] & ~(3 << (arrShift))) | (mask << (arrShift));
		}
	}
	shouldUpdatePackedGrid = false;
	return packedGrid_;
}

void Grid::fill(const GridCell cell) {
	std::fill(grid.get(), grid.get() + size, cell);
	shouldUpdatePackedGrid = true;
}

GridCell Grid::cellAtIndex(const uint32_t index) const {
	return grid.get()[index];
}

void Grid::setCellAtIndex(const uint32_t index, GridCell cell) {
	grid.get()[index] = cell;
	shouldUpdatePackedGrid = true;
}

void Grid::updateGeneration() {
	if (tasks == NULL) {
		tasks.reset(new std::unique_ptr<Task<TaskData>>[numberOfTasks]);
		const uint32_t taskGridHeight = height / numberOfTasks;
		for (uint32_t i = 1; i < numberOfTasks - 1; i++) {
			const uint32_t length = width * taskGridHeight;
			const uint32_t startRow = taskGridHeight * i;
			const auto output = new GridCell[length]{ GridCell::DEAD };
			const auto gridLocal = new GridCell[length + width * 2]{ GridCell::DEAD };
			tasks.get()[i] = std::unique_ptr<Task<TaskData>>(
				new Task<TaskData>{
					threadUpdateGrid,
					TaskData{
						*this,
						std::unique_ptr<GridCell[]>(gridLocal),
						grid,
						gridBuffer,
						std::unique_ptr<GridCell[]>(output),
						startRow,
						taskGridHeight,
						interrupt_flag
					}
				}
			);
		}
		{//top
			const uint32_t length = width * taskGridHeight;
			const auto output = new GridCell[length]{ GridCell::DEAD };
			const auto gridLocal = new GridCell[length + width * 2]{ GridCell::DEAD };
			tasks.get()[0] = std::unique_ptr<Task<TaskData>>(
				new Task<TaskData>{
					threadUpdateGridTop,
					TaskData{
						*this,
						std::unique_ptr<GridCell[]>(gridLocal),
						grid,
						gridBuffer,
						std::unique_ptr<GridCell[]>(output),
						0,
						taskGridHeight,
						interrupt_flag
					}
				}
			);
		}
		{ //bottom
			const uint32_t totalRows = (numberOfTasks - 1) * taskGridHeight;
			const uint32_t remainingRows = height - totalRows;
			const uint32_t length = width * remainingRows;
			const uint32_t startRows = height - remainingRows;
			const auto output = new GridCell[length]{ GridCell::DEAD };
			const auto gridLocal = new GridCell[length + width * 2]{ GridCell::DEAD };
			tasks.get()[numberOfTasks - 1] = std::unique_ptr<Task<TaskData>>(
				new Task<TaskData>{
					threadUpdateGridBot,
					TaskData{
						*this,
						std::unique_ptr<GridCell[]>(gridLocal),
						grid,
						gridBuffer,
						std::unique_ptr<GridCell[]>(output),
						startRows,
						remainingRows,
						interrupt_flag
					}
				}
			);
		}

		for (uint32_t i = 0; i < numberOfTasks; i++) {
			tasks.get()[i]->start();
		}
	}
	for (uint32_t i = 0; i < numberOfTasks; i++) {
		tasks.get()[i]->waitForResult();
	}
	for (uint32_t i = 0; i < numberOfTasks; i++) {
		const auto &data = tasks.get()[i]->data;
		for (int32_t i = 0; i < data.rowCount * width; i++) {
			grid[i + data.startRow * width] = data.output[i];
		}
	}

	for (uint32_t i = 1; i < numberOfTasks-1; i++) {
		const auto& data = tasks.get()[i]->data;
		const auto &gridLocal = data.gridLocal;
		const auto startRow = data.startRow;
		const auto rowCount = data.rowCount;
		std::memcpy(gridLocal.get(), &(grid.get()[(startRow-1)*width]), sizeof(GridCell) * (rowCount+2) * width);
		tasks.get()[i]->start();
	}
	{//top
		const auto& data = tasks.get()[0]->data;
		const auto& gridLocal = data.gridLocal;
		const auto rowCount = data.rowCount;
		std::memcpy(gridLocal.get(), &(grid.get()[size - width]), sizeof(GridCell) * 1 * width);
		std::memcpy(gridLocal.get(), &(grid.get()[0]), sizeof(GridCell) * (rowCount + 1) * width);
		tasks.get()[0]->start();
	}
	{//bottom
		const uint32_t i = numberOfTasks - 1;
		const auto& data = tasks.get()[i]->data;
		const auto& gridLocal = data.gridLocal;
		const auto startRow = data.startRow;
		const auto rowCount = data.rowCount;
		std::memcpy(gridLocal.get(), &(grid.get()[0]), sizeof(GridCell) * 1 * width);
		std::memcpy(gridLocal.get(), &(grid.get()[(startRow - 1) * width]), sizeof(GridCell) * (rowCount + 1) * width);
		tasks.get()[i]->start();
	}
	shouldUpdatePackedGrid = true;
}

GridCell Grid::updatedCell(const int32_t index) const {
	GridCell cell = grid[index];

	if (cell != GridCell::WALL) {
		uint32_t aliveNeighbours = 0;

		for (int yo = -1; yo <= 1; yo++) {
			for (int xo = -1; xo <= 1; xo++) {
				if (xo != 0 || yo != 0) {
					int x = ((index + xo) + width) % width,
						y = (((index / width) + yo) + height) % height;
					int offsetedIndex = x + width * y;
					if (cellAtIndex(offsetedIndex) == GridCell::ALIVE)
						aliveNeighbours++;
				}
			}
		}

		if (cell == GridCell::DEAD && aliveNeighbours == 3) {
			return GridCell::ALIVE;
		}
		else if (cell == GridCell::ALIVE && (aliveNeighbours < 2 || aliveNeighbours > 3)) {
			return GridCell::DEAD;
		}
	}

	return cell;
}

//void Grid::updateGeneration() {
//	for (int32_t i = 0; i < size; i++) {
//		gridBuffer[i] = updatedCell(i);
//	}
//
//	grid.swap(gridBuffer);
//	shouldUpdatePackedGrid = true; //mutable
//}

/*void Grid::updateGeneration() {
	const auto updatedCell = [](GridCell cell, unsigned char aliveNeighbours) {
		if (cell == GridCell::DEAD && aliveNeighbours == 3) {
			return GridCell::ALIVE;
		}
		else if (cell == GridCell::ALIVE && (aliveNeighbours < 2 || aliveNeighbours > 3)) {
			return GridCell::DEAD;
		}
		return cell;
	};

	const auto isCell = [this](size_t index) -> bool { return this->grid.get()[index] == GridCell::ALIVE; };
	const int lastElement = width - 1;
	const int lastRow = width * (height - 1);

	for (int row = 1; row < height - 1; row++) {
		const int rowIndex = width * row;
		for (int col = 1; col < width - 1; col++) {
			const int index = col + rowIndex;
			if (grid[index] != GridCell::WALL) {
				const unsigned char bottomRowNeighbours = isCell(index + width - 1) + isCell(index + width) + isCell(index + width + 1);
				const unsigned char topRowNeighbours =    isCell(index - width - 1) + isCell(index - width) + isCell(index - width + 1);
				const unsigned char aliveNeighbours = topRowNeighbours + isCell(index - 1) + isCell(index + 1) + bottomRowNeighbours;

				gridBuffer[index] = updatedCell(grid[index], aliveNeighbours);
			}
		}
		//first element
		{ 
			const int index = rowIndex;
			if (grid[index] != GridCell::WALL) {
				const unsigned char bottomRowNeighbours = isCell(index + width + lastElement) + isCell(index + width) + isCell(index + width + 1);
				const unsigned char topRowNeighbours =    isCell(index - width + lastElement) + isCell(index - width) + isCell(index - width + 1);
				const unsigned char aliveNeighbours = topRowNeighbours + isCell(index + lastElement) + isCell(index + 1) + bottomRowNeighbours;

				gridBuffer[rowIndex] = updatedCell(grid[rowIndex], aliveNeighbours);
			}
		}
		{ //last element
			const int index = rowIndex + lastElement;
			if (grid[index] != GridCell::WALL) {
				const unsigned char bottomRowNeighbours = isCell(index + width - 1) + isCell(index + width) + isCell(rowIndex + width);
				const unsigned char    topRowNeighbours = isCell(index - width - 1) + isCell(index - width) + isCell(rowIndex - width);
				const unsigned char aliveNeighbours = topRowNeighbours + isCell(index - 1) + isCell(rowIndex + 0) + bottomRowNeighbours;

				gridBuffer[rowIndex + lastElement] = updatedCell(grid[rowIndex + lastElement], aliveNeighbours);
			}
		}
	}

	{//first row
		for (int col = 1; col < width - 1; col++) {
			int index = col;
			if (grid[index] != GridCell::WALL) {
				const unsigned char topRowNeighbours = isCell(lastRow + index - 1) + isCell(lastRow + index) + isCell(lastRow + index + 1);
				const unsigned char bottomRowNeighbours = isCell(width + index - 1) + isCell(index + width) + isCell(width + index + 1);
				const unsigned char aliveNeighbours = topRowNeighbours + isCell(index - 1) + isCell(index + 1) + bottomRowNeighbours;

				gridBuffer[index] = updatedCell(grid[index], aliveNeighbours);
			}
		}

		{ //first element
			if (grid[0] != GridCell::WALL) {
				const unsigned char    topRowNeighbours = isCell(lastRow + lastElement) + isCell(lastRow + 0) + isCell(lastRow + 1);
				const unsigned char bottomRowNeighbours = isCell(width + lastElement)   + isCell(width)       + isCell(width + 1);
				const unsigned char aliveNeighbours = topRowNeighbours + isCell(lastElement) + isCell(1) + bottomRowNeighbours;

				gridBuffer[0] = updatedCell(grid[0], aliveNeighbours);
			}
		}
		{ //last element
			int index = lastElement;
			if (grid[index] != GridCell::WALL) {
				const unsigned char bottomRowNeighbours = isCell(lastRow + lastElement - 1) + isCell(lastRow + lastElement) + isCell(lastRow);
				const unsigned char    topRowNeighbours = isCell(width + lastElement - 1) + isCell(width + lastElement) + isCell(width);
				const unsigned char aliveNeighbours = topRowNeighbours + isCell(lastElement - 1) + isCell(0) + bottomRowNeighbours;

				gridBuffer[index] = updatedCell(grid[index], aliveNeighbours);
			}
		}
	}
	{//last row
		for (int col = 1; col < width - 1; col++) {
			int index = lastRow + col;
			if (grid[index] != GridCell::WALL) {
				const unsigned char topRowNeighbours = isCell(index - width - 1) + isCell(index - width) + isCell(index - width + 1);
				const unsigned char bottomRowNeighbours = isCell(col - 1) + isCell(col) + isCell(col + 1);
				const unsigned char aliveNeighbours = topRowNeighbours + isCell(lastRow + col - 1) + isCell(lastRow + col + 1) + bottomRowNeighbours;

				gridBuffer[index] = updatedCell(grid[index], aliveNeighbours);
			}
		}

		{ //first element
			int index = lastRow;
			if (grid[index] != GridCell::WALL) {
				const unsigned char    topRowNeighbours = isCell(index - width + lastElement) + isCell(index - width) + isCell(index - width + 1);
				const unsigned char bottomRowNeighbours = isCell(lastElement) + isCell(0) + isCell(1);
				const unsigned char aliveNeighbours = topRowNeighbours + isCell(index + lastElement) + isCell(index + 1) + bottomRowNeighbours;

				gridBuffer[index] = updatedCell(grid[index], aliveNeighbours);
			}
		}
		{ //last element
			int index = lastRow + lastElement;
			if (grid[index] != GridCell::WALL) {
				const unsigned char    topRowNeighbours = isCell(index - width - 1) + isCell(index - width) + isCell(lastRow - width);
				const unsigned char bottomRowNeighbours = isCell(lastElement - 1) + isCell(lastElement) + isCell(0);
				const unsigned char aliveNeighbours = topRowNeighbours + isCell(index - 1) + isCell(lastRow) + bottomRowNeighbours;

				gridBuffer[index] = updatedCell(grid[index], aliveNeighbours);
			}
		}
	}

	grid.swap(gridBuffer);
	shouldUpdatePackedGrid = true;
}*/
