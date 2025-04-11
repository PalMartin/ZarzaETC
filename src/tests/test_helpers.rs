#[macro_export]
macro_rules! assert_f64_eq {
    ($a:expr, $b:expr) => {
        if ($a - $b).abs() >= 1e-13 {
            panic!(
                "Assertion failed: `({} ≈ {})` (difference: {}, allowed error: 1e-13)",
                $a,
                $b,
                ($a - $b).abs()
            );
        }
    };
}

#[macro_export]
macro_rules! assert_f32_eq {
    ($a:expr, $b:expr) => {
        if ($a - $b).abs() >= 1e-6 {
            panic!(
                "Assertion failed: `({} ≈ {})` (difference: {}, allowed error: 1e-6)",
                $a,
                $b,
                ($a - $b).abs()
            );
        }
    };
}

