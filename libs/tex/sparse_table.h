/*
 * Copyright (C) 2015, Nils Moehrle
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#ifndef TEX_SPARSETABLE_HEADER
#define TEX_SPARSETABLE_HEADER

#include <vector>

#include <fstream>
#include <cstring>
#include <cerrno>
#include <iostream>

#include "util/file_system.h"
#include "util/exception.h"

#define HEADER "SPT"
#define VERSION "0.2"

/**
  * Class representing a sparse table optimized for row and column wise access.
  */
template <typename C, typename R, typename T>
class SparseTable {
    public:
        typedef std::vector<std::pair<R, T> > Column;
        typedef std::vector<std::pair<C, T> > Row;

    private:
        std::vector<Column> column_wise_data;
        std::vector<Row> row_wise_data;

        std::size_t nnz;
    public:
        SparseTable();
        SparseTable(C cols, R rows);

        C cols() const;
        R rows() const;

        Column const & col(C id) const;
        Row const & row(R id) const;

        void set_value(C col, R row, T value);

        std::size_t get_nnz(void) const;

        /**
          * Saves the SparseTable to the file given by filename.
          * Version and size information is stored in ascii, values in binary.
          * @throws util::FileException
          */
        static void save_to_file(SparseTable const & sparse_table, std::string const & filename);

        /**
          * Loads a SparseTable from the file given by filename.
          * @throws util::FileException if the file does not exist or if the header does not matches.
          */
        static void load_from_file(std::string const & filename, SparseTable * sparse_table);
};

template <typename C, typename R, typename T> std::size_t
SparseTable<C, R, T>::get_nnz(void) const {
    return nnz;
}

template <typename C, typename R, typename T> C
SparseTable<C, R, T>::cols() const {
    return column_wise_data.size();
}

template <typename C, typename R, typename T> R
SparseTable<C, R, T>::rows() const {
    return row_wise_data.size();
}

template <typename C, typename R, typename T> typename SparseTable<C, R, T>::Column const &
SparseTable<C, R, T>::col(C id) const {
    return column_wise_data[id];
}

template <typename C, typename R, typename T> typename SparseTable<C, R, T>::Row const &
SparseTable<C, R, T>::row(R id) const {
    return row_wise_data[id];
}

template <typename C, typename R, typename T>
SparseTable<C, R, T>::SparseTable() {
    nnz = 0;
}

template <typename C, typename R, typename T>
SparseTable<C, R, T>::SparseTable(C cols, R rows) {
    column_wise_data.resize(cols);
    row_wise_data.resize(rows);
    nnz = 0;
}

template <typename C, typename R, typename T> void
SparseTable<C, R, T>::set_value(C col, R row, T value) {
    column_wise_data[col].push_back(std::pair<R, T>(row, value));
    row_wise_data[row].push_back(std::pair<C, T>(col, value));
    nnz++;
}

template <typename C, typename R, typename T> void
SparseTable<C, R, T>::save_to_file(SparseTable const & sparse_table, const std::string &filename) {
    std::ofstream out(filename.c_str());
    if (!out.good())
        throw util::FileException(filename, std::strerror(errno));

    C cols = sparse_table.cols();
    R rows = sparse_table.rows();
    std::size_t nnz = sparse_table.get_nnz();
    out << HEADER << " " << VERSION << " " << cols << " " << rows << " " << nnz << std::endl;

    for (C col = 0; col < cols; ++col) {
        SparseTable::Column column = sparse_table.col(col);
        for (std::size_t i = 0; i < column.size(); ++i) {
            std::pair<R, T> entry = column[i];
            R row = entry.first;
            T value = entry.second;

            out.write((char const*)&col, sizeof(C));
            out.write((char const*)&row, sizeof(R));
            out.write((char const*)&value, sizeof(T));
        }
    }
    out.close();
}

template <typename C, typename R, typename T> void
SparseTable<C, R, T>::load_from_file(const std::string & filename, SparseTable<C, R, T> * sparse_table) {
    std::ifstream in(filename.c_str());
    if (!in.good())
        throw util::FileException(filename, std::strerror(errno));

    std::string header;

    in >> header;

    if (header != HEADER) {
        in.close();
        throw util::FileException(filename, "Not a SparseTable file!");
    }

    std::string version;
    in >> version;

    if (version != VERSION) {
        in.close();
        throw util::FileException(filename, "Incompatible version of SparseTable file!");
    }

    C cols;
    R rows;
    std::size_t nnz;
    in >> cols >> rows >> nnz;

    if (cols != sparse_table->cols() || rows != sparse_table->rows()) {
        in.close();
        throw util::FileException(filename, "SparseTable has different dimension!");
    }

    std::string buffer;

    /* Discard the rest of the line. */
    std::getline(in, buffer);

    C col;
    R row;
    T value;
    for (std::size_t i = 0; i < nnz; ++i) {
        in.read((char*)&col, sizeof(C));
        in.read((char*)&row, sizeof(R));
        in.read((char*)&value, sizeof(T));
        sparse_table->set_value(col, row, value);
    }

    in.close();
}

#endif /* TEX_SPARSETABLE_HEADER */
