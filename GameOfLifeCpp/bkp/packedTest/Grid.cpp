#include "Misc.h"
#include "Grid.h"

struct PackedIndex {
	uint32_t index, shift;

	PackedIndex(uint32_t index_, uint32_t shift_) : index(index_), shift(shift_) {}
};

static uint32_t cellMask = 0b11;
static PackedIndex packIndex(uint32_t index) {
	const uint32_t arrIndex = index / (32 / 2);
	const uint32_t arrShift = (index % (32 / 2)) * 2;
	return PackedIndex(arrIndex, arrShift);
}

static void setGridCellAtIndex(const std::shared_ptr<uint32_t[]> &array, const uint32_t index, GridCell cell) {
	const auto pi = packIndex(index);
	array.get()[pi.index] = array.get()[pi.index] & (~(cellMask << pi.shift)) | (uint32_t(cell) << pi.shift);
}

GridCell Grid::updatedCell(const uint32_t index) const {
	GridCell cell = cellAtIndex(index);

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

void Grid::updateGeneration() {
	for (int32_t i = 0; i < size; i++) {
		setGridCellAtIndex(gridBuffer, i, updatedCell(i));
	}

	grid.swap(gridBuffer);
	shouldUpdatePackedGrid = true;
}

const std::shared_ptr<const uint32_t[]> Grid::packedGrid() {
	return grid;
}

/*
void Grid::updateGeneration() {
	const auto updatedCell = [](GridCell cell, unsigned char aliveNeighbours) {
		if (cell == GridCell::DEAD && aliveNeighbours == 3) {
			return GridCell::ALIVE;
		}
		else if (cell == GridCell::ALIVE && (aliveNeighbours < 2 || aliveNeighbours > 3)) {
			return GridCell::DEAD;
		}
		return cell;
	};
	const auto gr = (grid.get());
	const auto isCell = [&gr](size_t index) -> bool { return gr[index] == GridCell::ALIVE; };
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

		{ //first element
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

void Grid::fill(const GridCell cell) {
	uint32_t mask = 0;
	for (unsigned char i = 0; i < 16; i++) {
		mask |= uint32_t(cell) << (i * 2);
	}
	std::fill(grid.get(), grid.get() + packedGridSize, mask);
	shouldUpdatePackedGrid = true;
}

GridCell Grid::cellAtIndex(const uint32_t index) const {
	const auto pi = packIndex(index);
	return static_cast<GridCell>(((grid.get()[pi.index] >> pi.shift) & cellMask));
}

void Grid::setCellAtIndex(const uint32_t index, GridCell cell) {
	setGridCellAtIndex(grid, index, cell);
	shouldUpdatePackedGrid = true;
}



