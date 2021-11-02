#pragma once

#include<iostream>
#include <stdint.h>
#include<memory>
#include "Vector.h"
#include <thread>
#include <atomic>
#include "Task.h"
#include <atomic>
#include <vector>
#include"MedianCounter.h"

enum class FieldCell : uint8_t
{
	DEAD = 0,
	ALIVE = 1,
	WALL = 0b10000
};

namespace fieldCell {
	bool isAlive(const FieldCell cell);
	bool isDead(const FieldCell cell);
	bool isWall(const FieldCell cell);
	FieldCell nextGeneration(const FieldCell cell, const uint32_t aliveNegihboursCount);
	inline constexpr char const* const asString(const FieldCell cell) {
		switch (cell) {
			case FieldCell::DEAD : return "dead" ;
			case FieldCell::ALIVE: return "alive";
			case FieldCell::WALL : return "wall" ;
			default: return "error";
		}
	}
}

class Field final {
public:
	class FieldPimpl;
	struct GridData;
	struct PackedGridData;
//public:
//	const uint32_t packedGridSize;
private:
	std::unique_ptr<FieldPimpl> gridPimpl;
	bool shouldUpdateGrid;
	bool isStopped;
	bool isFieldGPUBufferWriteOffset;
	std::atomic_bool gpuBufferLock_flag;

	const uint32_t numberOfTasks;
	std::unique_ptr<std::unique_ptr<Task<std::unique_ptr<GridData>>>[/*numberOfTasks*/]> gridTasks;
	std::atomic_bool interrupt_flag;
	std::vector<uint32_t> indecesToBrokenCells;

	const GLuint bufferP;
public:
	Field(const uint32_t gridWidth, const uint32_t gridHeight, const size_t numberOfTasks_, const GLuint bufferP, GLFWwindow *window, bool deployTasks);
	~Field();

	Field(Field const&) = delete;
	Field& operator=(Field const&) = delete;
public:
	void updateGeneration();
	void fill(const FieldCell cell);

	FieldCell cellAtIndex(const uint32_t index) const;

	void setCellAtIndex(const uint32_t index, FieldCell cell);
	void setCellAtCoord(const vec2i& coord, FieldCell cell);

	FieldCell cellAtCoord(const vec2i& coord) const;
	FieldCell cellAtCoord(const int32_t column, const int32_t row) const;

	vec2i indexAsCoord(const int32_t index) const;

	uint32_t coordAsIndex(const vec2i& coord) const;
	uint32_t coordAsIndex(const int32_t column, const int32_t row) const;
	
	uint32_t normalizeIndex(const int32_t index) const;
	vec2i normalizeCoord(const vec2i& coord) const;

	//delete &, append const
	bool &isGridUpdated();

	uint32_t width() const;
	uint32_t height() const;
	uint32_t size() const;

	FieldCell* grid();

	void stopAllGridTasks();
	void startAllGridTasks();

	bool isFieldBufferWriteOffset();
private:
	void waitForGridTasks();
	void deployGridTasks();
	uint32_t currentWriteOffset();
	uint32_t currentReadOffset();
};

inline bool Field::isFieldBufferWriteOffset() {
	return isFieldGPUBufferWriteOffset;
}


inline void Field::setCellAtCoord(const vec2i& coord, FieldCell cell) {
	setCellAtIndex(coordAsIndex(coord), cell);
}

inline FieldCell Field::cellAtCoord(const vec2i& coord) const {
	return cellAtIndex(static_cast<int32_t>(coordAsIndex(coord)));
}

inline FieldCell Field::cellAtCoord(const int32_t column, const int32_t row) const {
	return cellAtCoord(vec2i(column, row));
}

inline vec2i Field::indexAsCoord(const int32_t index) const {
	const auto index_n = normalizeIndex(index);
	const auto x = index_n % width();
	const auto y = index_n / width();
	return vec2i(x, y);
}

inline uint32_t Field::coordAsIndex(const vec2i& coord) const {
	const auto coord_n = normalizeCoord(coord);
	return coord_n.x + coord_n.y * width();
}

inline uint32_t Field::coordAsIndex(const int32_t column, const int32_t row) const {
	return coordAsIndex(vec2i(column, row));
}

inline uint32_t Field::normalizeIndex(const int32_t index) const {
	return misc::umod(index, static_cast<uint32_t>(size()));
}
inline vec2i Field::normalizeCoord(const vec2i& coord) const {
	return vec2i(misc::mod(coord.x, width()), misc::mod(coord.y, height()));
}

inline bool &Field::isGridUpdated() {
	return shouldUpdateGrid;
}