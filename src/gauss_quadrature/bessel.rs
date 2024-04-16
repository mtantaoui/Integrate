use std::f64::consts::PI;

use num::{one, Float, ToPrimitive, Unsigned};

/// Computes the kth zero of the $J_0(x)$ Bessel function.
///
/// # Notes
///
/// Note that the first 20 zeros are tabulated.  After that, they are computed
fn bessel_j0<F: Float, U: Unsigned + ToPrimitive>(k: U) -> f64 {
    const J_Z: [f64; 20] = [
        2.40482555769577276862163187933E+00,
        5.52007811028631064959660411281E+00,
        8.65372791291101221695419871266E+00,
        11.7915344390142816137430449119E+00,
        14.9309177084877859477625939974E+00,
        18.0710639679109225431478829756E+00,
        21.2116366298792589590783933505E+00,
        24.3524715307493027370579447632E+00,
        27.4934791320402547958772882346E+00,
        30.6346064684319751175495789269E+00,
        33.7758202135735686842385463467E+00,
        36.9170983536640439797694930633E+00,
        40.0584257646282392947993073740E+00,
        43.1997917131767303575240727287E+00,
        46.3411883716618140186857888791E+00,
        49.4826098973978171736027615332E+00,
        52.6240518411149960292512853804E+00,
        55.7655107550199793116834927735E+00,
        58.9069839260809421328344066346E+00,
        62.0484691902271698828525002646E+00,
    ];

    let r: f64;
    let mut r2: f64;
    let mut z: f64;

    let mut tmp: f64;

    if J_Z.len() > 20 {
        z = PI * (k.to_f64().unwrap() - 0.25E+00);
        r = 1.0E+00 / z;
        r2 = r * r;

        tmp = r2 * 0.509225462402226769498681286758E+08;
        tmp += -0.849353580299148769921876983660E+06;
        tmp *= r2;
        tmp += 0.186904765282320653831636345064E+05;
        tmp *= r2;
        tmp += -0.567644412135183381139802038240E+03;
        tmp *= r2;
        tmp += 0.253364147973439050099206349206E+02;
        tmp *= r2;
        tmp += -0.182443876720610119047619047619E+01;
        tmp *= r2;
        tmp += 0.246028645833333333333333333333E+00;
        tmp *= r2;
        tmp += 0.125E+00;
        tmp *= r;
        z += tmp;
    } else {
        z = J_Z[k.to_usize().unwrap() - 1];
    }

    z
}

fn formula(x: f64, x2: f64) -> f64 {
    x * (0.202642367284675542887758926420E+00
        + x2 * x2
            * (-0.303380429711290253026202643516E-03
                + x2 * (0.198924364245969295201137972743E-03
                    + x2 * (-0.228969902772111653038747229723E-03
                        + x2 * (0.433710719130746277915572905025E-03
                            + x2 * (-0.123632349727175414724737657367E-02
                                + x2 * (0.496101423268883102872271417616E-02
                                    + x2 * (-0.266837393702323757700998557826E-01
                                        + x2 * (0.185395398206345628711318848386E+00)))))))))
}

/// Computes the kth zero of the $J_0(x)$ Bessel function.
///
/// # Notes
///
/// Note that the first 20 zeros are tabulated.  After that, they are computed
fn bessel_j1_squared<F: Float, U: Unsigned + ToPrimitive>(k: U) -> f64 {
    const J_1: &[f64; 21] = &[
        0.269514123941916926139021992911E+00,
        0.115780138582203695807812836182E+00,
        0.0736863511364082151406476811985E+00,
        0.0540375731981162820417749182758E+00,
        0.0426614290172430912655106063495E+00,
        0.0352421034909961013587473033648E+00,
        0.0300210701030546726750888157688E+00,
        0.0261473914953080885904584675399E+00,
        0.0231591218246913922652676382178E+00,
        0.0207838291222678576039808057297E+00,
        0.0188504506693176678161056800214E+00,
        0.0172461575696650082995240053542E+00,
        0.0158935181059235978027065594287E+00,
        0.0147376260964721895895742982592E+00,
        0.0137384651453871179182880484134E+00,
        0.0128661817376151328791406637228E+00,
        0.0120980515486267975471075438497E+00,
        0.0114164712244916085168627222986E+00,
        0.0108075927911802040115547286830E+00,
        0.0102603729262807628110423992790E+00,
        0.00976589713979105054059846736696E+00,
    ];

    let x: f64;
    let x2: f64;
    let z: f64;

    let mut tmp: f64;

    if J_1.len() < k.to_usize().unwrap() {
        x = 1.0 / (k.to_f64().unwrap() - 0.25);
        x2 = x * x;

        z = formula(x, x2);
    } else {
        z = J_1[k.to_usize().unwrap() - 1];
    }
    z
}

/// Computes the $K^{th}$ pair of an $N$-point Gauss-Legendre rule.
///
/// # Discussion
///
/// $\theta$ values of the zeros are in $\[0,pi\]$, and monotonically increasing.
///
fn glpair<U: Unsigned>(n: U, k: U) {}

/// Computes the $K^{th}$ pair of an $N$-point Gauss-Legendre rule.
///
/// # Discussion
///
/// $\theta$ values of the zeros are in $\[0,pi\]$, and monotonically increasing.
///
pub fn glpairs<U: Unsigned + ToPrimitive + PartialOrd + Copy>(n: U, k: U) {
    if n < one::<U>() {
        panic!("GLPAIRS - FATAL ERROR \n Illegal value of N");
    }

    if k < one::<U>() || n < k {
        panic!("GLPAIRS - FATAL ERROR \n Illegal value of K");
    }

    let kcopy;
    if n < k + k - one::<U>() {
        kcopy = n - k + one();
    } else {
        kcopy = k;
    }

    // get the bessel zero
    let w = 1.0 / (n.to_f64().unwrap() + 0.5);
    let nu = bessel_j0::<f64, U>(kcopy);
    let theta = w * nu;
    let y = theta * theta;

    // get the asymptotic BesselJ(1, nu) squared
    let b = bessel_j1_squared::<f64, U>(kcopy);

    // get chebyshev interpolants for nodes

    let sf1t = (((((-1.29052996274280508473467968379E-12 * y
        + 2.40724685864330121825976175184E-10)
        * y
        - 3.13148654635992041468855740012E-08)
        * y
        + 0.275573168962061235623801563453E-05)
        * y
        - 0.148809523713909147898955880165E-03)
        * y
        + 0.416666666665193394525296923981E-02)
        * y
        - 0.416666666666662959639712457549E-01;
}
