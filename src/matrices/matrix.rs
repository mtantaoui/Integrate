use core::fmt;
use std::cmp::min;
use std::ops::Add;
use std::{fmt::Debug, marker::Send};

use num::{Float, Zero};
use rayon::prelude::*;

#[derive(Clone)]
pub struct FloatMatrix<F: Float> {
    data: Vec<F>, // matrix elements stored in Row-major Order
    nrows: usize, // number of rows
    ncols: usize, // number of columns
}

impl<F: Float + Sized + Send + Debug> FloatMatrix<F> {
    pub fn new(data: Vec<F>, nrows: usize, ncols: usize) -> FloatMatrix<F> {
        if data.len() != nrows * ncols {
            panic!("creation failed");
        }
        FloatMatrix { data, nrows, ncols }
    }

    fn get_element(&self, i: usize, j: usize) -> F {
        let index = i * self.ncols + j;
        let msg = format!("Failed to get element ({},{})", i, j);
        *self.data.get(index).expect(&msg)
    }

    pub fn set_element(&mut self, i: usize, j: usize, new_element: F) {
        self.data[i * self.ncols + j] = new_element;
    }

    fn size(&self) -> usize {
        self.nrows * self.ncols
    }

    pub fn nrows(&self) -> usize {
        self.nrows
    }

    pub fn ncols(&self) -> usize {
        self.ncols
    }

    pub fn zero(nrows: usize, ncols: usize) -> FloatMatrix<F> {
        let size = nrows * ncols;
        let mut data: Vec<F> = Vec::new();
        (0..size)
            .into_par_iter()
            .map(|_| F::zero())
            .collect_into_vec(&mut data);
        FloatMatrix { nrows, ncols, data }
    }

    pub fn transpose(&mut self) {
        let (nrows, ncols) = (self.nrows, self.ncols);
        (self.nrows, self.ncols) = (ncols, nrows);

        for i in 0..self.nrows {
            for j in 0..min(i, self.ncols) {
                let previous_index = i * self.ncols + j;
                let new_index = j * self.nrows + i;

                self.data.swap(previous_index, new_index);
            }
        }
    }
}

impl<F: Float + Sized + Send + Sync + Debug> Zero for FloatMatrix<F> {
    fn zero() -> Self {
        unimplemented!()
    }

    fn is_zero(&self) -> bool {
        is_zero(self.data.as_ref())
    }

    fn set_zero(&mut self) {
        self.data = vec![F::zero(); self.size()]
    }
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
        }
    }
}

impl<F: Float + fmt::Debug + Send> fmt::Display for FloatMatrix<F> {
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
