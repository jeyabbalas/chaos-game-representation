use chaos_game_representation::Point;


const EPSILON: f64 = 0.0000001;

#[test]
pub fn add_test() {
    let sum = Point { x: 0.5, y: 1.3 } + Point { x: 2.0, y: 3.7 };
    assert!(((sum.x - 2.5) as f64).abs() < EPSILON);
    assert!(((sum.y - 5.0) as f64).abs() < EPSILON);
}


#[test]
pub fn subtract_test() {
    let sub = Point { x: 0.5, y: 3.7 } - Point { x: 2.0, y: 1.3 };
    assert!(((sub.x + 1.5) as f64).abs() < EPSILON);
    assert!(((sub.y - 2.4) as f64).abs() < EPSILON);
}


#[test]
pub fn divide_test() {
    let div = Point { x: 2.0, y: 0.5 } / 2.0;
    assert!(((div.x - 1.0) as f64).abs() < EPSILON);
    assert!(((div.y - 0.25) as f64).abs() < EPSILON);
}