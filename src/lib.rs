#![feature(generic_const_exprs)]
#![allow(incomplete_features)]

extern crate nalgebra as na;

pub mod builder;
mod kernel;
mod powers;
pub mod rbf;
mod test;
