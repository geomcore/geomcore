//! [`Continuity`]: how smooth a curve or a surface direction is at a
//! parameter, as its knots guarantee it.

/// Continuity of a curve or surface direction at a parameter, as its
/// knots guarantee it.
///
/// `Order(k)` means the shape is at least `C^k` there: its derivatives up
/// to order `k` agree on both sides. A B-spline knot reports
/// `Order(degree − multiplicity)`. This is the order
/// the knots guarantee, which the shape may exceed (a knot inserted into a
/// smooth curve lowers the reported order, not the curve's smoothness).
/// `Smooth` means no knot constrains the shape there: every derivative
/// agrees on both sides.
///
/// Every `Order(k)` sorts below `Smooth`, and orders compare by `k`, so
/// [`Ord::min`] picks the weaker of two continuities.
///
/// # Examples
///
/// ```
/// use geomcore::Continuity;
///
/// assert!(Continuity::Order(1) < Continuity::Order(2));
/// assert!(Continuity::Order(u32::MAX) < Continuity::Smooth);
/// assert_eq!(Continuity::Order(2).min(Continuity::Smooth), Continuity::Order(2));
/// assert!(Continuity::Smooth.is_at_least(2));
/// assert!(!Continuity::Order(1).is_at_least(2));
/// ```
// The derived ordering relies on `Order` being declared before `Smooth`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Continuity {
    /// At least `C^k`: derivatives up to order `k` agree on both sides.
    Order(u32),
    /// No knot constrains the shape: every derivative agrees on both
    /// sides.
    Smooth,
}

impl Continuity {
    /// The finite order `k` of `Order(k)`, or `None` when smooth.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::Continuity;
    ///
    /// assert_eq!(Continuity::Order(2).order(), Some(2));
    /// assert_eq!(Continuity::Smooth.order(), None);
    /// ```
    #[must_use]
    #[inline]
    pub const fn order(self) -> Option<u32> {
        match self {
            Continuity::Order(order) => Some(order),
            Continuity::Smooth => None,
        }
    }

    /// Whether the shape is at least `C^order` here: `Smooth`, or
    /// `Order(k)` with `k >= order`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::Continuity;
    ///
    /// assert!(Continuity::Order(2).is_at_least(2));
    /// assert!(!Continuity::Order(0).is_at_least(1));
    /// assert!(Continuity::Smooth.is_at_least(u32::MAX));
    /// ```
    #[must_use]
    #[inline]
    pub const fn is_at_least(self, order: u32) -> bool {
        match self {
            Continuity::Order(k) => k >= order,
            Continuity::Smooth => true,
        }
    }
}
