extern crate terrain_erosion; // Replace 'my_library' with the name of your library crate

use terrain_erosion::adjust_outflow_to_total_water; // Import library functions as needed
use terrain_erosion::calculate_k_factor; // Import library functions as needed
use terrain_erosion::calculate_terrain_gradient; // Import library functions as needed
use terrain_erosion::calculate_water_outflow_flux; // Import library functions as needed
use terrain_erosion::scale_flux_with_k_factor; // Import library functions as needed
use terrain_erosion::update_bottomward_flux; // Import library functions as needed
use terrain_erosion::update_leftward_flux; // Import library functions as needed
use terrain_erosion::update_rightward_flux; // Import library functions as needed
use terrain_erosion::update_topward_flux; // Import library functions as needed
use terrain_erosion::CellProperties; // Import library functions as needed

fn main() {
    let mut cell = CellProperties {
        terrain_height: 3.14 * 2.0,
        water_height: 5.0,
        suspended_sediment: 2.0,
        water_outflow_flux: (0.0, 0.0, 0.0, 0.0), // Sample values for f_L, f_R, f_T, f_B
        velocity: (0.0, 0.0),                     // Sample velocity
    };
    let surrounding_cell = CellProperties {
        terrain_height: 0.0,
        water_height: 5.0,
        suspended_sediment: 2.0,
        water_outflow_flux: (0.0, 0.0, 0.0, 0.0), // Sample values for f_L, f_R, f_T, f_B
        velocity: (0.0, 0.0),                     // Sample velocity
    };
    let left_cell = surrounding_cell.clone();
    let right_cell = surrounding_cell.clone();
    let top_cell = surrounding_cell.clone();
    let bottom_cell = surrounding_cell.clone();
    let terrain_gradient = calculate_terrain_gradient(
        &cell,
        &left_cell,
        &right_cell,
        &top_cell,
        &bottom_cell,
    );

    // Update the leftward water outflow flux
    let delta_time = 1.0;
    let area = 1.0;
    let length = 1.0;
    let gravity = 9.81;
    let delta_height_left = cell.terrain_height - left_cell.terrain_height;
    let delta_height_right = cell.terrain_height - right_cell.terrain_height;
    let delta_height_top = cell.terrain_height - top_cell.terrain_height;
    let delta_height_bottom = cell.terrain_height - bottom_cell.terrain_height;
    let k_factor = 0.8;

    let water_outflow_flux = calculate_water_outflow_flux(
        &cell,
        &terrain_gradient,
        5.0,
        delta_time,
        area,
        gravity,
        k_factor,
    );
    cell.water_outflow_flux = water_outflow_flux;

    println!("Water Outflow Flux: {:?}", water_outflow_flux);

    let new_leftward_flux = update_leftward_flux(
        &cell,
        delta_time,
        area,
        length,
        gravity,
        delta_height_left,
        k_factor,
    );

    let new_rightward_flux = update_rightward_flux(
        &cell,
        delta_time,
        area,
        length,
        gravity,
        delta_height_right,
        k_factor,
    );

    let new_topward_flux = update_topward_flux(
        &cell,
        delta_time,
        area,
        length,
        gravity,
        delta_height_top,
        k_factor,
    );

    let new_bottomward_flux = update_bottomward_flux(
        &cell,
        delta_time,
        area,
        length,
        gravity,
        delta_height_bottom,
        k_factor,
    );

    // Calculate K factor
    let k_factor = calculate_k_factor(
        cell.water_height,
        delta_time,
        new_leftward_flux,
        new_rightward_flux,
        new_topward_flux,
        new_bottomward_flux,
    );

    // Apply K factor to all flux components
    let scaled_leftward_flux = scale_flux_with_k_factor(new_leftward_flux, k_factor);

    let scaled_rightward_flux = scale_flux_with_k_factor(new_rightward_flux, k_factor);

    let scaled_topward_flux = scale_flux_with_k_factor(new_topward_flux, k_factor);

    let scaled_bottomward_flux = scale_flux_with_k_factor(new_bottomward_flux, k_factor);

    // Check and adjust the total outflow to not exceed the total water content
    adjust_outflow_to_total_water(
        &mut cell,
        scaled_leftward_flux,
        scaled_rightward_flux,
        scaled_topward_flux,
        scaled_bottomward_flux,
    );

    println!("Updated Cell: {:#?}", cell);
}
