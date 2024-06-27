use core::fmt;
use itertools::concat;
use rayon::iter::{IndexedParallelIterator, IntoParallelIterator, ParallelIterator};
use std::ops::Add;
use std::{fmt::Debug, marker::Send};

use num::Float;

#[derive(Clone)]
pub enum MatrixStorageType {
    RowMajorOrder,
    ColumnMajorOrder,
}

pub trait Matrix<F: Float> {
    fn new(data: Vec<F>, nrows: usize, ncols: usize, storage_type: MatrixStorageType) -> Self;
    fn get_element(&self, i: usize, j: usize) -> F;
    fn get_index(&self, i: usize, j: usize) -> usize;
    fn set_element(&mut self, i: usize, j: usize, new_element: F);
    fn size(&self) -> usize;
    fn nrows(&self) -> usize;
    fn ncols(&self) -> usize;
    fn zero(nrows: usize, ncols: usize, storage_type: MatrixStorageType) -> Self;
    fn transpose(&mut self);
    fn is_zero(&self) -> bool;
    fn set_zero(&mut self);
}

#[derive(Clone)]
pub struct FloatMatrix<F: Float> {
    data: Vec<F>, // matrix elements stored in Row-major Order
    nrows: usize, // number of rows
    ncols: usize, // number of columns
    storage_type: MatrixStorageType,
}

impl<F: Float + Sized + Send + Debug + Sync> Matrix<F> for FloatMatrix<F> {
    fn new(
        data: Vec<F>,
        nrows: usize,
        ncols: usize,
        storage_type: MatrixStorageType,
    ) -> FloatMatrix<F> {
        if data.len() != nrows * ncols {
            panic!("creation failed");
        }
        FloatMatrix {
            data,
            nrows,
            ncols,
            storage_type,
        }
    }

    fn get_element(&self, i: usize, j: usize) -> F {
        let index = match self.storage_type {
            MatrixStorageType::RowMajorOrder => i * self.ncols + j,
            MatrixStorageType::ColumnMajorOrder => j * self.nrows + i,
        };

        let msg = format!("Failed to get element ({},{})", i, j);
        *self.data.get(index).expect(&msg)
    }

    fn get_index(&self, i: usize, j: usize) -> usize {
        match self.storage_type {
            MatrixStorageType::RowMajorOrder => i * self.ncols + j,
            MatrixStorageType::ColumnMajorOrder => j * self.nrows + i,
        }
    }

    fn set_element(&mut self, i: usize, j: usize, new_element: F) {
        let index = self.get_index(i, j);
        self.data[index] = new_element;
    }

    fn size(&self) -> usize {
        self.nrows * self.ncols
    }

    fn nrows(&self) -> usize {
        self.nrows
    }

    fn ncols(&self) -> usize {
        self.ncols
    }

    fn zero(nrows: usize, ncols: usize, storage_type: MatrixStorageType) -> FloatMatrix<F> {
        let size = nrows * ncols;
        let mut data: Vec<F> = Vec::new();
        (0..size)
            .into_par_iter()
            .map(|_| F::zero())
            .collect_into_vec(&mut data);
        FloatMatrix {
            nrows,
            ncols,
            data,
            storage_type,
        }
    }

    // works for row major column order only !!!
    fn transpose(&mut self) {
        self.data = match self.storage_type {
            MatrixStorageType::RowMajorOrder => {
                transpose_row_major_order(self.data.as_ref(), self.nrows, self.ncols)
            }
            MatrixStorageType::ColumnMajorOrder => {
                transpose_column_major_order(self.data.as_ref(), self.nrows, self.ncols)
            }
        };

        let (nrows, ncols) = (self.nrows, self.ncols);
        (self.nrows, self.ncols) = (ncols, nrows);
    }

    fn is_zero(&self) -> bool {
        is_zero(self.data.as_ref())
    }

    fn set_zero(&mut self) {
        self.data = vec![F::zero(); self.size()]
    }
    //     // // works only for square matrices
    //     // pub fn transpose(&mut self) {
    //     //     let (nrows, ncols) = (self.nrows, self.ncols);
    //     //     (self.nrows, self.ncols) = (ncols, nrows);

    //     //     for i in 0..self.nrows {
    //     //         for j in 0..min(i, self.ncols) {
    //     //             let previous_index = i * self.ncols + j;
    //     //             let new_index = j * self.nrows + i;

    //     //             self.data.swap(previous_index, new_index);
    //     //         }
    //     //     }
    //     // }
}

impl<F: Float + Sized + Sync + Send> Add for FloatMatrix<F> {
    type Output = FloatMatrix<F>;

    fn add(self, rhs: Self) -> Self::Output {
        if self.nrows != rhs.nrows || self.ncols != rhs.ncols {
            panic!("failed")
        }

        let data1: &[F] = self.data.as_ref();
        let data2: &[F] = rhs.data.as_ref();

        let data = add(data1, data2);

        FloatMatrix {
            data,
            nrows: self.nrows,
            ncols: self.ncols,
            storage_type: MatrixStorageType::RowMajorOrder,
        }
    }
}

impl<F: Float + fmt::Debug + Send + Sync> fmt::Display for FloatMatrix<F> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut output = "".to_owned();
        for i in 0..self.nrows {
            for j in 0..self.ncols {
                let element_output = format!("{:?}\t", self.get_element(i, j));
                output.push_str(element_output.as_ref());
            }
            output.push('\n');
        }
        write!(f, "{}", output)
    }
}

fn add<F: Float + Send + Sync>(data1: &[F], data2: &[F]) -> Vec<F> {
    if data1.len() != data2.len() {
        panic!("crashed")
    }
    let size = data1.len();

    (0..size)
        .into_par_iter()
        .zip_eq(data1)
        .zip_eq(data2)
        .map(|((_, &element1), &element2)| element1 + element2)
        .collect()
}

fn is_zero<F: Float + Send + Sync>(data: &[F]) -> bool {
    let size = data.len();

    (0..size)
        .into_par_iter()
        .zip(data)
        .all(|(_, &element)| element.is_zero())
}

fn transpose_row_major_order<F: Float + Send + Sync>(
    data: &[F],
    nrows: usize,
    ncols: usize,
) -> Vec<F> {
    (0..ncols)
        .into_par_iter()
        .map(|j| {
            let mut row = Vec::new();
            for i in 0..nrows {
                let element = data[i * ncols + j];
                row.push(element);
            }
            row
        })
        .reduce(|| Vec::new(), |acc, row| concat(vec![acc, row]))
}

fn transpose_column_major_order<F: Float + Send + Sync>(
    data: &[F],
    nrows: usize,
    ncols: usize,
) -> Vec<F> {
    (0..nrows)
        .into_par_iter()
        .map(|j| {
            let mut column = Vec::new();

            for i in 0..ncols {
                let element = data[i + j * nrows];
                column.push(element);
            }

            column
        })
        .reduce(|| Vec::new(), |acc, column| concat(vec![acc, column]))
}
