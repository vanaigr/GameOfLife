#include "Misc.h"
#include "Grid.h"

static GridCell updatedCell(const Grid& grid, const int32_t index) {
	GridCell cell = grid.cellAtIndex(index);
	//const vec2i cellCoord = grid.indexAsCoord(index);

	if (cell != GridCell::WALL) {
		uint32_t aliveNeighbours = 0;

		for (int yo = -1; yo <= 1; yo++) {
			for (int xo = -1; xo <= 1; xo++) {
				if (xo != 0 || yo != 0) {
					int x = ((index + xo) + grid.width) % grid.width,
						y = (((index / grid.width) + yo) + grid.height) % grid.height;
					int offsetedIndex = x + grid.width * y;
					//const auto offset = vec2i{ xo, yo };
					//const auto coord = cellCoord + offset;
					//if (grid.cellAtCoord(coord) == GridCell::ALIVE)
					if (grid.cellAtIndex(offsetedIndex) == GridCell::ALIVE)
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

void Grid::updateGeneration() {
	for (int32_t i = 0; i < size; i++) {
		gridBuffer[i] = updatedCell(*this, i);
	}

	grid.swap(gridBuffer);
	shouldUpdatePackedGrid = true; //mutable
}

const std::shared_ptr<const uint32_t[]> Grid::packedGrid() const {
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
	shouldUpdatePackedGrid = false; //mutable
	return packedGrid_;
}


//void Grid::updateGeneration() {
//	const auto updatedCell = [](GridCell cell, unsigned char aliveNeighbours) {
//		if (cell == GridCell::DEAD && aliveNeighbours == 3) {
//			return GridCell::ALIVE;
//		}
//		else if (cell == GridCell::ALIVE && (aliveNeighbours < 2 || aliveNeighbours > 3)) {
//			return GridCell::DEAD;
//		}
//		return cell;
//	};
//	const auto gr = (grid.get());
//	const auto isCell = [&gr](size_t index) -> bool { return gr[index] == GridCell::ALIVE; };
//	const int lastElement = width - 1;
//	const int lastRow = width * (height - 1);
//
//	for (int row = 1; row < height - 1; row++) {
//		const int rowIndex = width * row;
//		for (int col = 1; col < width - 1; col++) {
//			const int index = col + rowIndex;
//			if (grid[index] != GridCell::WALL) {
//				const unsigned char bottomRowNeighbours = isCell(index + width - 1) + isCell(index + width) + isCell(index + width + 1);
//				const unsigned char topRowNeighbours =    isCell(index - width - 1) + isCell(index - width) + isCell(index - width + 1);
//				const unsigned char aliveNeighbours = topRowNeighbours + isCell(index - 1) + isCell(index + 1) + bottomRowNeighbours;
//
//				gridBuffer[index] = updatedCell(grid[index], aliveNeighbours);
//			}
//		}
//
//		{ //first element
//			const int index = rowIndex;
//			if (grid[index] != GridCell::WALL) {
//				const unsigned char bottomRowNeighbours = isCell(index + width + lastElement) + isCell(index + width) + isCell(index + width + 1);
//				const unsigned char topRowNeighbours =    isCell(index - width + lastElement) + isCell(index - width) + isCell(index - width + 1);
//				const unsigned char aliveNeighbours = topRowNeighbours + isCell(index + lastElement) + isCell(index + 1) + bottomRowNeighbours;
//
//				gridBuffer[rowIndex] = updatedCell(grid[rowIndex], aliveNeighbours);
//			}
//		}
//		{ //last element
//			const int index = rowIndex + lastElement;
//			if (grid[index] != GridCell::WALL) {
//				const unsigned char bottomRowNeighbours = isCell(index + width - 1) + isCell(index + width) + isCell(rowIndex + width);
//				const unsigned char    topRowNeighbours = isCell(index - width - 1) + isCell(index - width) + isCell(rowIndex - width);
//				const unsigned char aliveNeighbours = topRowNeighbours + isCell(index - 1) + isCell(rowIndex + 0) + bottomRowNeighbours;
//
//				gridBuffer[rowIndex + lastElement] = updatedCell(grid[rowIndex + lastElement], aliveNeighbours);
//			}
//		}
//	}
//	{//first row
//		for (int col = 1; col < width - 1; col++) {
//			int index = col;
//			if (grid[index] != GridCell::WALL) {
//				const unsigned char topRowNeighbours = isCell(lastRow + index - 1) + isCell(lastRow + index) + isCell(lastRow + index + 1);
//				const unsigned char bottomRowNeighbours = isCell(width + index - 1) + isCell(index + width) + isCell(width + index + 1);
//				const unsigned char aliveNeighbours = topRowNeighbours + isCell(index - 1) + isCell(index + 1) + bottomRowNeighbours;
//
//				gridBuffer[index] = updatedCell(grid[index], aliveNeighbours);
//			}
//		}
//
//		{ //first element
//			if (grid[0] != GridCell::WALL) {
//				const unsigned char    topRowNeighbours = isCell(lastRow + lastElement) + isCell(lastRow + 0) + isCell(lastRow + 1);
//				const unsigned char bottomRowNeighbours = isCell(width + lastElement)   + isCell(width)       + isCell(width + 1);
//				const unsigned char aliveNeighbours = topRowNeighbours + isCell(lastElement) + isCell(1) + bottomRowNeighbours;
//
//				gridBuffer[0] = updatedCell(grid[0], aliveNeighbours);
//			}
//		}
//		{ //last element
//			int index = lastElement;
//			if (grid[index] != GridCell::WALL) {
//				const unsigned char bottomRowNeighbours = isCell(lastRow + lastElement - 1) + isCell(lastRow + lastElement) + isCell(lastRow);
//				const unsigned char    topRowNeighbours = isCell(width + lastElement - 1) + isCell(width + lastElement) + isCell(width);
//				const unsigned char aliveNeighbours = topRowNeighbours + isCell(lastElement - 1) + isCell(0) + bottomRowNeighbours;
//
//				gridBuffer[index] = updatedCell(grid[index], aliveNeighbours);
//			}
//		}
//	}
//	{//last row
//		for (int col = 1; col < width - 1; col++) {
//			int index = lastRow + col;
//			if (grid[index] != GridCell::WALL) {
//				const unsigned char topRowNeighbours = isCell(index - width - 1) + isCell(index - width) + isCell(index - width + 1);
//				const unsigned char bottomRowNeighbours = isCell(col - 1) + isCell(col) + isCell(col + 1);
//				const unsigned char aliveNeighbours = topRowNeighbours + isCell(lastRow + col - 1) + isCell(lastRow + col + 1) + bottomRowNeighbours;
//
//				gridBuffer[index] = updatedCell(grid[index], aliveNeighbours);
//			}
//		}
//
//		{ //first element
//			int index = lastRow;
//			if (grid[index] != GridCell::WALL) {
//				const unsigned char    topRowNeighbours = isCell(index - width + lastElement) + isCell(index - width) + isCell(index - width + 1);
//				const unsigned char bottomRowNeighbours = isCell(lastElement) + isCell(0) + isCell(1);
//				const unsigned char aliveNeighbours = topRowNeighbours + isCell(index + lastElement) + isCell(index + 1) + bottomRowNeighbours;
//
//				gridBuffer[index] = updatedCell(grid[index], aliveNeighbours);
//			}
//		}
//		{ //last element
//			int index = lastRow + lastElement;
//			if (grid[index] != GridCell::WALL) {
//				const unsigned char    topRowNeighbours = isCell(index - width - 1) + isCell(index - width) + isCell(lastRow - width);
//				const unsigned char bottomRowNeighbours = isCell(lastElement - 1) + isCell(lastElement) + isCell(0);
//				const unsigned char aliveNeighbours = topRowNeighbours + isCell(index - 1) + isCell(lastRow) + bottomRowNeighbours;
//
//				gridBuffer[index] = updatedCell(grid[index], aliveNeighbours);
//			}
//		}
//	}
//
//	grid.swap(gridBuffer);
//}

void Grid::fill(const GridCell cell) {
	std::fill(grid.get(), grid.get() + size, cell);
	shouldUpdatePackedGrid = true; //mutable
}

GridCell Grid::cellAtIndex(const int32_t index) const {
	auto index_n = normalizeIndex(index);
	return grid.get()[index_n];
}

void Grid::setCellAtIndex(const int32_t index, GridCell cell) {
	auto index_n = normalizeIndex(index);
	grid.get()[index_n] = cell;
	shouldUpdatePackedGrid = true; //mutable
}