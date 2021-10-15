#pragma once

#include<iostream>
#include <stdint.h>
#include<memory>
#include "Vector.h"
#include <thread>
#include <atomic>
#include "Task.h"
#include <atomic>

enum class GridCell : unsigned char
{
	DEAD,
	ALIVE,
	WALL
};

class Grid;

struct TaskData {
	Grid& gridO;
	std::unique_ptr<GridCell[]> output;
	uint32_t startRow;
	uint32_t rowCount;
	std::unique_ptr < GridCell[]> gridLocal; //size of output + 2 rows, index + row
//private:
	std::unique_ptr < GridCell[]>& grid;
	std::unique_ptr < GridCell[]>& gridBuffer;//
	std::atomic_bool &interrupt_flag;//

public:
	TaskData(
		Grid& gridO_,
		std::unique_ptr < GridCell[]> gridLocal_,
		std::unique_ptr < GridCell[]>& grid_,
		std::unique_ptr < GridCell[]>& gridBuffer_,
		std::unique_ptr<GridCell[]>&& output_,
		uint32_t startRow_,
		uint32_t rowCount_,
		std::atomic_bool& interrupt_flag_

	) :
		gridO(gridO_),
		gridLocal(std::move(gridLocal_)),
		grid(grid_),
		gridBuffer(gridBuffer_),
		output(std::move(output_)),
		startRow(startRow_),
		rowCount(rowCount_),
		interrupt_flag(interrupt_flag_)
	{}

	const std::unique_ptr<GridCell[]>& currentGrid(bool isBuffering) {
		if (isBuffering) return gridBuffer;
		return grid;
	}

	bool shouldInterrupt(bool isBuffering) {
		return isBuffering && interrupt_flag.load();
	}
};

class Grid final {
public:
	const uint32_t width;
	const uint32_t height;
	const uint32_t size;
	const uint32_t packedGridSize;
private:
	std::unique_ptr<GridCell[/*size*/]> grid;
	std::unique_ptr<GridCell[/*size*/]> gridBuffer;
	const std::shared_ptr<uint32_t[/*packedGridSize*/]> packedGrid_;
	bool shouldUpdatePackedGrid;

	const uint32_t numberOfTasks;
	std::unique_ptr<std::unique_ptr<Task<TaskData>>[/*numberOfThreads*/]> tasks;
	std::atomic_bool interrupt_flag;

public:
	Grid(const uint32_t gridWidth, const uint32_t gridHeight, const size_t numberOfTasks_);
	~Grid() = default;

	Grid(Grid const&) = delete;
	Grid& operator=(Grid const&) = delete;
public:
	void updateGeneration();
	void fill(const GridCell cell);

	GridCell cellAtIndex(const uint32_t index) const;

	void setCellAtIndex(const uint32_t index, GridCell cell);
	void setCellAtCoord(const vec2i& coord, GridCell cell);

	GridCell cellAtCoord(const vec2i& coord) const;
	GridCell cellAtCoord(const int32_t column, const int32_t row) const;

	vec2i indexAsCoord(const int32_t index) const;

	int32_t coordAsIndex(const vec2i& coord) const;
	int32_t coordAsIndex(const int32_t column, const int32_t row) const;
	
	int32_t normalizeIndex(const int32_t index) const;
	vec2i normalizeCoord(const vec2i& coord) const;

	const std::shared_ptr<const uint32_t[]> packedGrid();
	bool isPackedGridUpdated() const;
private:
	GridCell updatedCell(const int32_t index) const;
};

inline void Grid::setCellAtCoord(const vec2i& coord, GridCell cell) {
	setCellAtIndex(coordAsIndex(coord), cell);
}

inline GridCell Grid::cellAtCoord(const vec2i& coord) const {
	return cellAtIndex(static_cast<int32_t>(coordAsIndex(coord)));
}

inline GridCell Grid::cellAtCoord(const int32_t column, const int32_t row) const {
	return cellAtCoord(vec2i(column, row));
}

inline vec2i Grid::indexAsCoord(const int32_t index) const {
	const auto index_n = normalizeIndex(index);
	const auto x = index_n % width;
	const auto y = index_n / width;
	return vec2i(x, y);
}

inline int32_t Grid::coordAsIndex(const vec2i& coord) const {
	const auto coord_n = normalizeCoord(coord);
	return coord_n.x + coord_n.y * width;
}

inline int32_t Grid::coordAsIndex(const int32_t column, const int32_t row) const {
	return coordAsIndex(vec2i(column, row));
}

inline int32_t Grid::normalizeIndex(const int32_t index) const {
	return misc::mod(index, size);
}
inline vec2i Grid::normalizeCoord(const vec2i& coord) const {
	return vec2i(misc::mod(coord.x, width), misc::mod(coord.y, height));
}

inline bool Grid::isPackedGridUpdated() const {
	return shouldUpdatePackedGrid;
}