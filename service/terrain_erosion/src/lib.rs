mod terrain_erosion {
    extern crate rand;

    use rand::Rng;

    pub fn simulate_rainfall(heightmap: &mut Vec<Vec<f64>>, x: usize, y: usize, intensity: f64) {
        let erosion_depth = 0.1;
        let mut current_x = x;
        let mut current_y = y;
    
        let mut remaining_intensity = intensity;
    
        while remaining_intensity > 0.0 {
            // Calculate the steepest downhill neighbor
            let neighbor = find_downhill_neighbor(heightmap, current_x, current_y);
            
            // // Calculate the height difference between the current and next position
            // let height_difference = heightmap[current_x][current_y] - heightmap[next_x][next_y];
    
            // // Calculate the amount of erosion based on the height difference
            // let erosion_depth = calculate_erosion_depth(height_difference, remaining_intensity);
    
            // // Update the heightmap to simulate erosion
            // heightmap[current_x][current_y] -= erosion_depth;
    
            // // Move to the next position
            // current_x = next_x;
            // current_y = next_y;
    
            // // Reduce the remaining rainfall intensity
            remaining_intensity -= erosion_depth;
        }
    }

    pub fn find_downhill_neighbor(
        height_map: &Vec<Vec<f64>>,
        x: usize,
        y: usize,
    ) -> Option<(usize, usize)> {
        let width = height_map.len();
        let height = height_map[0].len();
    
        let mut lowest_neighbor = None;
        let mut lowest_elevation = f64::MAX;
    
        for dx in -1..=1 {
            for dy in -1..=1 {
                // Calculate neighbor coordinates
                let nx = x as isize + dx;
                let ny = y as isize + dy;
                // Check if the neighbor is within bounds
                if nx >= 0 && ny >= 0 && nx < width as isize && ny < height as isize {

                    // Get the elevation of the neighbor
                    let neighbor_elevation = height_map[nx as usize][ny as usize];
    
                    // Compare with the lowest elevation found so far
                    if neighbor_elevation < lowest_elevation {
                        lowest_elevation = neighbor_elevation;
                        lowest_neighbor = Some((nx as usize, ny as usize));
                    }
                }
            }
        }
    
        lowest_neighbor
    }

    // Define other erosion functions and types
}

pub use terrain_erosion::*; // Re-export the module's contents

#[cfg(test)]
mod tests {
    use super::terrain_erosion::find_downhill_neighbor;
    use super::terrain_erosion::simulate_rainfall;

    #[test]
    fn test_simulate_rainfall() {
        // Arrange
        let mut height_map: Vec<Vec<f64>> = vec![
            vec![1.0, 2.0, 3.0],
            vec![4.0, 0.0, 5.0],
            vec![6.0, 7.0, 8.0],
        ];
        let initial_height_map = height_map.clone();
        // Act
        simulate_rainfall(&mut height_map, 0, 0, 1.0);

        // Assert
        assert_eq!(initial_height_map, height_map);
    }

    #[test]
    fn test_find_downhill_neighbor() {
        // Arrange
        let height_map: Vec<Vec<f64>> = vec![
            vec![1.0, 2.0, 3.0],
            vec![4.0, 0.0, 5.0],
            vec![6.0, 7.0, 8.0],
        ];

        // Act/Assert
        assert_eq!(find_downhill_neighbor(&height_map, 0, 0), Some((1, 1)));
        assert_eq!(find_downhill_neighbor(&height_map, 0, 1), Some((1, 1)));
        assert_eq!(find_downhill_neighbor(&height_map, 0, 2), Some((1, 1)));

        assert_eq!(find_downhill_neighbor(&height_map, 1, 0), Some((1, 1)));
        assert_eq!(find_downhill_neighbor(&height_map, 1, 1), Some((1, 1)));
        assert_eq!(find_downhill_neighbor(&height_map, 1, 2), Some((1, 1)));

        assert_eq!(find_downhill_neighbor(&height_map, 2, 0), Some((1, 1)));
        assert_eq!(find_downhill_neighbor(&height_map, 2, 1), Some((1, 1)));
        assert_eq!(find_downhill_neighbor(&height_map, 2, 2), Some((1, 1)));
    }
}
