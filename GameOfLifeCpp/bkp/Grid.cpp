#include "Misc.h"
#include "Grid.h"

GridCell updatedCell(const Grid& grid, const uint32_t index) {
	GridCell cell = grid.cellAtIndex(index);
	const vec2i cellCoord = grid.indexAsCoord(index);

	if (cell != GridCell::WALL) {
		uint32_t aliveNeighbours = 0;

		for (int yo = -1; yo <= 1; yo++) {
			for (int xo = -1; xo <= 1; xo++) {
				if (xo != 0 || yo != 0) {
					//int x = ((index + xo) + width) % gridWidth,
					//	y = (((index / gridWidth) + yo) + gridHeight) % gridHeight;
					const auto offset = vec2i{ xo, yo };
					const auto coord = cellCoord + offset;
					if (grid.cellAtCoord(coord) == GridCell::ALIVE)
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
	for (size_t i = 0; (int32_t)i < size; i++) {
		gridBuffer[i] = updatedCell(*this, i);
	}

	grid.swap(gridBuffer);
}

void Grid::fill(const GridCell cell) {
	std::fill(grid.get(), grid.get() + size, GridCell::DEAD);
}

GridCell Grid::cellAtIndex(const int32_t index) const {
	auto index_n = normalizeIndex(index);
	return grid[index_n];
}

void Grid::setCellAtIndex(const int32_t index, GridCell cell) {
	auto index_n = normalizeIndex(index);
	grid[index_n] = cell;
}