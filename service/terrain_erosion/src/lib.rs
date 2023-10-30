pub mod properties;

pub mod terrain_erosion {
    #![allow(unused_mut, unused_variables, unused_imports)]
    extern crate rand;

    use super::properties::CellProperties;
    use rand::Rng;

    pub fn simulate_rainfall(
        heightmap: &Vec<Vec<f64>>,
        cell_properties: &mut Vec<Vec<CellProperties>>,
        rain_rate: f64,
        time_step: f64,
        kr_parameter: f64,
    ) {
        println!("Simulating Rainfall");
        let rows = heightmap.len();
        let cols = heightmap[0].len();

        for x in 0..rows {
            for y in 0..cols {
                let dt = cell_properties[x][y].water_height;
                let d1 = dt + time_step * rain_rate * kr_parameter;
                // Update the water_height property of the corresponding cell
                cell_properties[x][y].water_height = d1;
            }
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

    // fₗ₊∆ₜᵗ = max(0, fₗᵗ(x, y) + ∆ₜ · A · (g · ∆hᴸ(x, y)/l))
    pub fn get_flux(
        current_properties_direction_flux: &f64,
        delta_time: &f64,
        area: &f64,
        length: &f64,
        gravity: &f64,
        height_difference: &f64,
    ) -> f64 {
        let flux = if current_properties_direction_flux.is_finite() {
            current_properties_direction_flux
                + delta_time * area * (gravity * height_difference / length)
        } else {
            delta_time * area * (gravity * height_difference / length)
        };
        f64::max(0.0, flux)
    }

    // K = max(1, d1 · lx · ly(fL + fR + fT + fB)·∆t)
    ///    K factor is the max of 1 or
    ///    initial water height
    ///    times the distances between grid cells
    ///    divided by the sum of all flux
    ///    times the change in time
    pub fn calculate_k_factor(
        water_height: &f64,
        delta_time: &f64,
        fluxes: &crate::properties::DirectionalProperties,
    ) -> f64 {
        let denominator = (fluxes.left + fluxes.right + fluxes.top + fluxes.bottom) * delta_time;
        f64::max(1.0, water_height / denominator)
    }

    // fₜ₊∆ₜⁱ(x, y) = K · fⁱₜ₊∆ₜ ,i = L,R,T,B.
    pub fn scale_flux_with_k_factor(
        flux: &crate::properties::DirectionalProperties,
        k_factor: &f64,
    ) -> crate::properties::DirectionalProperties {
        crate::properties::DirectionalProperties {
            left: k_factor * flux.left,
            right: k_factor * flux.right,
            top: k_factor * flux.top,
            bottom: k_factor * flux.bottom,
        }
    }

    // ∆V(x,y) = ∆t · (∑ fᵢₙ −∑ fₒᵤₜ)
    // = ∆t ·(
    //    fₜ₊∆ₜᴿ(x−1, y) +
    //    fₜ₊∆ₜᵀ(x, y−1) +
    //    fₜ₊∆ₜᴸ(x+1, y) +
    //    fₜ₊∆ₜᴮ(x, y+1)
    //    − ∑ ᵢ₌L,T,R,B fₜ₊∆ₜⁱ(x, y) // scale_flux_factor_with_k_factor
    //  )
    pub fn get_water_height_change(
        delta_time: &f64,
        scaled_inflow_fluxes: &crate::properties::DirectionalProperties,
        scaled_outflow_fluxes: &crate::properties::DirectionalProperties,
    ) -> f64 {
        delta_time
            * (scaled_inflow_fluxes.left
                + scaled_inflow_fluxes.right
                + scaled_inflow_fluxes.top
                + scaled_inflow_fluxes.bottom
                - scaled_outflow_fluxes.left
                - scaled_outflow_fluxes.right
                - scaled_outflow_fluxes.top
                - scaled_outflow_fluxes.bottom)
    }

    #[derive(Debug)]
    pub struct TerrainHeightDifference {
        pub left: f64,
        pub right: f64,
        pub top: f64,
        pub bottom: f64,
    }

    pub fn calculate_height_differences(
        cell: &CellProperties,
        left_cell: &CellProperties,
        right_cell: &CellProperties,
        top_cell: &CellProperties,
        bottom_cell: &CellProperties,
    ) -> crate::properties::DirectionalProperties {
        let left_height_difference = (cell.terrain_height + cell.suspended_sediment)
            - (left_cell.terrain_height + left_cell.suspended_sediment);
        let right_height_difference = (cell.suspended_sediment + cell.terrain_height)
            - (right_cell.suspended_sediment + right_cell.terrain_height);
        let top_height_difference = (cell.terrain_height + cell.suspended_sediment)
            - (top_cell.suspended_sediment + top_cell.terrain_height);
        let bottom_height_difference = (cell.suspended_sediment + cell.terrain_height)
            - (bottom_cell.suspended_sediment + bottom_cell.terrain_height);

        crate::properties::DirectionalProperties {
            left: left_height_difference,
            right: right_height_difference,
            top: top_height_difference,
            bottom: bottom_height_difference,
        }
    }

    // α(x, y)
    /// Calculate the local tilt angle in radians.
    pub fn get_local_tilt_angle(
        current_cell: &CellProperties,
        leftward_cell: &CellProperties,
        rightward_cell: &CellProperties,
        topward_cell: &CellProperties,
        bottomward_cell: &CellProperties,
    ) -> f64 {
        let delta_h_x = rightward_cell.terrain_height - leftward_cell.terrain_height;
        let delta_h_y = topward_cell.terrain_height - bottomward_cell.terrain_height;

        let local_tilt_angle_rad = (delta_h_x.powi(2) + delta_h_y.powi(2)).sqrt().atan();

        local_tilt_angle_rad
    }

    // ∆Wₓ = ½(fᴿ(x-1, y) - fᴸ(x, y) + fᴿ(x, y) - fᴸ(x+1, y))
    pub fn get_x_hydraulic_erosion_and_deposition(
        scaled_outflow_fluxes: &crate::properties::DirectionalProperties,
        scaled_inflow_fluxes: &crate::properties::DirectionalProperties,
    ) -> f64 {
        let numerator: f64 = scaled_inflow_fluxes.left - scaled_outflow_fluxes.left
            + scaled_outflow_fluxes.right
            - scaled_inflow_fluxes.right;
        let denominator: f64 = 2.0;
        numerator / denominator
    }

    // ∆Wᵧ = ½(fᵀ(x, y-1) - fᴮ(x, y) + fᵀ(x, y) - fᴮ(x, y+1))
    pub fn get_y_hydraulic_erosion_and_deposition(
        scaled_outflow_fluxes: &crate::properties::DirectionalProperties,
        scaled_inflow_fluxes: &crate::properties::DirectionalProperties,
    ) -> f64 {
        let numerator: f64 = scaled_inflow_fluxes.top - scaled_outflow_fluxes.top
            + scaled_outflow_fluxes.bottom
            - scaled_inflow_fluxes.bottom;
        let denominator: f64 = 2.0;
        numerator / denominator
    }

    // C(x, y) = Kc · sin(α(x, y) · |−→v (x, y)| · lₘₐₓ)
    /// This function provides a value for sediment transfer.
    /// It takes into account a maximum sediment transfer at different water heights,
    /// but does not take into account 3D collision.
    pub fn get_transport_capacity(
        global_sediment_capacity: &f64,
        tilt_radians: &f64,
        h_erosion_and_deposition: &f64,
        v_erosion_and_deposition: &f64,
        optional_lmax: Option<&f64>,
    ) -> f64 {
        let lmax = optional_lmax.unwrap_or(&1.0);
        let water_velocity_magnitude = (h_erosion_and_deposition.powi(2) + v_erosion_and_deposition.powi(2)).sqrt();
        let sin_of = tilt_radians * water_velocity_magnitude;
        global_sediment_capacity * f64::sin(sin_of)
    }

    // lmax(x) =
    //  0, x ≤ 0
    //  1, x ≥ Kdmax
    //  1−(Kdmax −x)/Kdmax, 0 < x < Kdmax
    /// This function is used to provide a limiting factor to
    /// sediment transportation at different water heights.
    pub fn get_lmax(water_height: &f64, maximum_erosion_depth: &f64) -> f64 {
        if water_height <= &0.0 {
            0.0
        } else if water_height >= maximum_erosion_depth {
            1.0
        } else {
            1.0 - (maximum_erosion_depth - water_height) / maximum_erosion_depth
        }
    }

    // bt + ∆t = bt − ∆t · Rt(x, y) · Ks(C − st)
    // s1 = st + ∆t · Rt(x, y) · Ks(C − st)
    // d3 = d2 + ∆t · Rt(x, y) · Ks(C − st),
    //
    // bt + ∆t = bt + ∆t · Kd(st − C)
    // s1 = st − ∆t · Kd(st − C)
    // d3 = d2 − ∆t · Kd(st − C),
    /// This function will produce a tuple for new cell values.
    /// The original formula factors transportation rate into the 
    /// dissolve step. This can be addressed in future work, but is
    /// not required (Rₜ).
    pub fn get_updated_layers(
        cell: &CellProperties,
        capacity: &f64,
        dissolution_rate: &f64,
        desposition_rate: &f64,
        delta_time: &f64,
    ) -> (f64, f64, f64) {
        let sediment = &cell.suspended_sediment;
        let bed_thickness = &cell.terrain_height;
        let water_height = &cell.water_height;
        if sediment < capacity {
            // Dissolve some soil in water
            let dissolve_volume = delta_time * dissolution_rate * (capacity - sediment);
            let new_bed_thickness = bed_thickness - dissolve_volume;
            let new_sediment = sediment + dissolve_volume;
            let new_water_height = water_height + dissolve_volume;
    
            (new_sediment, new_bed_thickness, new_water_height)
        } else {
            // Dispose some of the transported sediment
            let dispose_volume = delta_time * desposition_rate * (sediment - capacity);
            let new_bed_thickness = bed_thickness + dispose_volume;
            let new_sediment = sediment - dispose_volume;
            let new_water_height = water_height - dispose_volume;
    
            (new_sediment, new_bed_thickness, new_water_height)
        }
    }
}

pub use terrain_erosion::*; // Re-export the module's contents

#[cfg(test)]
mod tests {
    use crate::properties::CellPropertiesBuilder;

    use super::terrain_erosion::adjust_outflow_to_total_water;

    #[test]
    fn test_adjust_outflow_to_total_water() {
        let (flux_l, flux_r, flux_t, flux_b) = (
            7.8480000000000008,
            7.8480000000000008,
            7.8480000000000008,
            7.8480000000000008,
        );
        let current_cell = CellPropertiesBuilder::new().water_height(1.0).build();
        assert_eq!(
            adjust_outflow_to_total_water(&current_cell, &flux_l, &flux_r, &flux_t, &flux_b),
            (1.0, 1.0, 1.0, 1.0)
        );
    }
}
