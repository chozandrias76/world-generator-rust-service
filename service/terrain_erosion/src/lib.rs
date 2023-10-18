mod terrain_erosion {
    #![allow(unused_mut, unused_variables, unused_imports)]
    extern crate rand;

    use rand::Rng;

    // Define a struct to represent cell properties
    #[derive(Debug, Clone)]
    pub struct CellProperties {
        pub terrain_height: f64,
        pub water_height: f64,
        pub suspended_sediment: f64,
        pub water_outflow_flux: (f64, f64, f64, f64), // (f_L, f_R, f_T, f_B)
        pub velocity: (f64, f64),
    }

    pub fn simulate_rainfall(heightmap: &mut Vec<Vec<f64>>, x: usize, y: usize, intensity: f64) {
        println!("Simulating Rainfall");
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
        println!(
            "Rain Fall preticipated from {intensity} to zero",
            intensity = intensity,
        )
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

    pub fn update_leftward_flux(
        current_properties: &CellProperties,
        delta_time: f64,
        area: f64,
        length: f64,
        gravity: f64,
        delta_height_left: f64,
        k_factor: f64,
    ) -> f64 {
        let current_leftward_flux = current_properties.water_outflow_flux.0;
        let new_leftward_flux = current_leftward_flux
            + delta_time * area * (gravity * delta_height_left * k_factor / length);
        f64::max(0.0, new_leftward_flux)
    }

    pub fn update_rightward_flux(
        current_properties: &CellProperties,
        delta_time: f64,
        area: f64,
        length: f64,
        gravity: f64,
        delta_height_right: f64,
        k_factor: f64,
    ) -> f64 {
        let current_rightward_flux = current_properties.water_outflow_flux.1;
        let new_rightward_flux = current_rightward_flux
            + delta_time * area * (gravity * delta_height_right * k_factor / length);
        f64::max(0.0, new_rightward_flux)
    }

    pub fn update_topward_flux(
        current_properties: &CellProperties,
        delta_time: f64,
        area: f64,
        length: f64,
        gravity: f64,
        delta_height_top: f64,
        k_factor: f64,
    ) -> f64 {
        let current_topward_flux = current_properties.water_outflow_flux.2;
        let new_topward_flux = current_topward_flux
            + delta_time * area * (gravity * delta_height_top * k_factor / length);
        f64::max(0.0, new_topward_flux)
    }

    pub fn update_bottomward_flux(
        current_properties: &CellProperties,
        delta_time: f64,
        area: f64,
        length: f64,
        gravity: f64,
        delta_height_bottom: f64,
        k_factor: f64,
    ) -> f64 {
        let current_bottomward_flux = current_properties.water_outflow_flux.3;
        let new_bottomward_flux = current_bottomward_flux
            + delta_time * area * (gravity * delta_height_bottom * k_factor / length);
        f64::max(0.0, new_bottomward_flux)
    }

    pub fn calculate_k_factor(
        water_height: f64,
        delta_time: f64,
        leftward_flux: f64,
        rightward_flux: f64,
        topward_flux: f64,
        bottomward_flux: f64,
    ) -> f64 {
        let denominator =
            (leftward_flux + rightward_flux + topward_flux + bottomward_flux) * delta_time;
        f64::max(1.0, water_height / denominator)
    }

    pub fn scale_flux_with_k_factor(flux: f64, k_factor: f64) -> f64 {
        flux * k_factor
    }

    pub fn adjust_outflow_to_total_water(
        current_properties: &mut CellProperties,
        leftward_flux: f64,
        rightward_flux: f64,
        topward_flux: f64,
        bottomward_flux: f64,
    ) {
        let total_flux = leftward_flux + rightward_flux + topward_flux + bottomward_flux;
        let total_water_content =
            current_properties.water_height + current_properties.suspended_sediment;
        println!("total_flux Cell: {}", total_flux);
        println!("total_water_content Cell: {}", total_water_content);

        if total_flux > total_water_content {
            let scaling_factor = total_water_content / total_flux;
            // Scale down all flux components to match the available water content
            current_properties.water_outflow_flux.0 *= scaling_factor;
            current_properties.water_outflow_flux.1 *= scaling_factor;
            current_properties.water_outflow_flux.2 *= scaling_factor;
            current_properties.water_outflow_flux.3 *= scaling_factor;
        }
    }

    pub struct TerrainGradient {
        left: f64,
        right: f64,
        top: f64,
        bottom: f64,
    }

    pub fn calculate_terrain_gradient(
        cell: &CellProperties,
        left_cell: &CellProperties,
        right_cell: &CellProperties,
        top_cell: &CellProperties,
        bottom_cell: &CellProperties,
    ) -> TerrainGradient {
        // Calculate terrain gradients in each direction
        let left_gradient = cell.terrain_height - left_cell.terrain_height;
        let right_gradient = right_cell.terrain_height - cell.terrain_height;
        let top_gradient = cell.terrain_height - top_cell.terrain_height;
        let bottom_gradient = bottom_cell.terrain_height - cell.terrain_height;

        // Calculate slopes (absolute values of gradients)
        let left_slope = left_gradient.abs();
        let right_slope = right_gradient.abs();
        let top_slope = top_gradient.abs();
        let bottom_slope = bottom_gradient.abs();

        TerrainGradient {
            left: left_slope,
            right: right_slope,
            top: top_slope,
            bottom: bottom_slope,
        }
    }

    pub fn calculate_water_outflow_flux(
        cell: &CellProperties,
        terrain_gradient: &TerrainGradient,
        water_height: f64,
        delta_time: f64,
        area: f64,
        gravity: f64,
        k_factor: f64,
    ) -> (f64, f64, f64, f64) {
        // Calculate the hydraulic radius (R) based on area and wetted perimeter
        let wetted_perimeter = 2.0 * (area + delta_time);
        let hydraulic_radius = area / wetted_perimeter;

        // Calculate the cross-sectional flow area (A)
        let cross_sectional_area = hydraulic_radius * wetted_perimeter;

        // Calculate the Manning's equation for flow rate
        let flow_rate = (1.0 / k_factor)
            * cross_sectional_area
            * (terrain_gradient.left.powf(5.0 / 3.0)
                + terrain_gradient.right.powf(5.0 / 3.0)
                + terrain_gradient.top.powf(5.0 / 3.0)
                + terrain_gradient.bottom.powf(5.0 / 3.0))
            .powf(3.0 / 5.0)
            * (1.0 / (1.0 - k_factor.powf(2.0)))
            * (water_height.powf(5.0 / 3.0) * gravity.powf(1.0 / 2.0) * terrain_gradient.left);

        // Distribute flow equally to all four directions
        let flux_left = flow_rate / 4.0;
        let flux_right = flow_rate / 4.0;
        let flux_top = flow_rate / 4.0;
        let flux_bottom = flow_rate / 4.0;

        (
            flux_left * delta_time,
            flux_right * delta_time,
            flux_top * delta_time,
            flux_bottom * delta_time,
        )
    }
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
