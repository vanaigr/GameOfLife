#pragma once

#include<iostream>
#include <stdint.h>
#include<memory>
#include "Vector.h"

enum GridCell : unsigned char
{
	DEAD,
	ALIVE,
	WALL
};

class Grid final {
public:
	const uint32_t width;
	const uint32_t height;
	const uint32_t size;
	const uint32_t packedGridSize;
private:
	std::shared_ptr<uint32_t[/*packedGridSize*/]> grid;
	std::shared_ptr<uint32_t[/*packedGridSize*/]> gridBuffer;
	bool shouldUpdatePackedGrid;
public:
	Grid(const uint32_t gridWidth, const uint32_t gridHeight) :
		width(gridWidth), height(gridHeight), size(gridWidth * gridHeight),
		packedGridSize((width* height) / (32 / 2) + 1),
		grid{ new uint32_t[packedGridSize]{ 0 } },
		gridBuffer{ new uint32_t[packedGridSize]{ 0 } },
		shouldUpdatePackedGrid(false) {}
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
	GridCell updatedCell(const uint32_t index) const;
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