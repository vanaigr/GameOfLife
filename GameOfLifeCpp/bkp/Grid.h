#pragma once

#include <stdint.h>
#include<memory>
#include "Vector.h"

enum class GridCell : unsigned char
{
	DEAD,
	ALIVE,
	WALL
};

class Grid final {
public:
	const int32_t width;
	const int32_t height;
	const int32_t size;
private:
	std::unique_ptr<GridCell[/*size*/]> grid;
	std::unique_ptr<GridCell[/*size*/]> gridBuffer;
public:
	Grid(const int32_t gridWidth, const int32_t gridHeight) : 
		width(gridWidth), height(gridHeight), size(gridWidth * gridHeight),
		grid { new GridCell[size]{ GridCell::DEAD } },
		gridBuffer{ new GridCell[size]{ GridCell::DEAD } } {}
	~Grid() = default;

	Grid(Grid const&) = delete;
	Grid& operator=(Grid const&) = delete;
public:
	void updateGeneration();
	void fill(const GridCell cell);

	GridCell cellAtIndex(const int32_t index) const;

	void setCellAtIndex(const int32_t index, GridCell cell);
	void setCellAtCoord(const vec2i& coord, GridCell cell);

	GridCell cellAtCoord(const vec2i& coord) const;
	GridCell cellAtCoord(const int32_t column, const int32_t row) const;

	vec2i indexAsCoord(const int32_t index) const;

	int32_t coordAsIndex(const vec2i& coord) const;
	int32_t coordAsIndex(const int32_t column, const int32_t row) const;
	
	int32_t normalizeIndex(const int32_t index) const;
	vec2i normalizeCoord(const vec2i& coord) const;
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