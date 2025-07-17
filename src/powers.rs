pub(crate) fn monomial_exponents<const VARIABLES: usize, const DEGREE: usize>()
-> Vec<[i32; VARIABLES]> {
    let mut result = Vec::new();

    fn rec<const VARIABLES: usize>(
        cur: &mut [i32; VARIABLES],
        remain: usize,
        var: usize,
        n: usize,
        out: &mut Vec<[i32; VARIABLES]>,
    ) {
        if var + 1 == n {
            cur[var] = remain as i32;
            out.push(cur.clone());
        } else {
            for i in 0..=remain {
                cur[var] = i as i32;
                rec(cur, remain - i, var + 1, n, out);
            }
        }
    }

    for total_deg in 0..=DEGREE {
        let mut cur = [0; VARIABLES];
        rec(&mut cur, total_deg, 0, VARIABLES, &mut result);
    }
    result
}

#[cfg(test)]
mod tests {
    use crate::powers::monomial_exponents;

    #[test]
    fn test_monomial_powers() {
        assert_eq!(monomial_exponents::<2, 0>(), vec![[0, 0]]);
        assert_eq!(monomial_exponents::<2, 1>(), vec![[0, 0], [0, 1], [1, 0]]);
        assert_eq!(
            monomial_exponents::<2, 2>(),
            vec![[0, 0], [0, 1], [1, 0], [0, 2], [1, 1], [2, 0]]
        );
        assert_eq!(
            monomial_exponents::<2, 3>(),
            vec![
                [0, 0],
                [0, 1],
                [1, 0],
                [0, 2],
                [1, 1],
                [2, 0],
                [0, 3],
                [1, 2],
                [2, 1],
                [3, 0]
            ]
        );
    }
}
