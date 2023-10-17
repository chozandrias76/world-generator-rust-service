use image::{ImageBuffer, Rgba};
use noise::utils::{NoiseMapBuilder, PlaneMapBuilder};
use noise::{Fbm, NoiseFn, Perlin, RidgedMulti};
use rand::Rng;

fn main() {
    let (map_width, map_height, seed, scale, octaves, persistence, lacunarity, offset) = (1000, 1000, 0, 200.0, 5, 0.5, 2.0, (0.0, 0.0));
    let noise_map = generate_noise_map(map_width, map_height, seed, scale, octaves, persistence, lacunarity, offset);
    let base_file_name = "output_heightmap";
    let output_file_name = format!("{}_{}_{}_{}_{}_{}_{}_{}_{}_{}.png",
        base_file_name, map_width, map_height, seed, scale, octaves, persistence, lacunarity, offset.0, offset.1
    );
    generate_colored_heightmap(&noise_map, &output_file_name);
}

pub fn interpolate_color(color1: Rgba<u8>, color2: Rgba<u8>, t: f64) -> Rgba<u8> {
    let r = ((1.0 - t) * color1.0[0] as f64 + t * color2.0[0] as f64) as u8;
    let g = ((1.0 - t) * color1.0[1] as f64 + t * color2.0[1] as f64) as u8;
    let b = ((1.0 - t) * color1.0[2] as f64 + t * color2.0[2] as f64) as u8;
    Rgba([r, g, b, 255])
}

pub fn generate_colored_heightmap(heightmap: &Vec<Vec<f64>>, output_path: &str) {
    let width = heightmap.len();
    let height = heightmap[0].len();

    let mut image_buffer = ImageBuffer::new(width as u32, height as u32);

    for x in 0..width {
        for y in 0..height {
            let height_value = heightmap[x][y];

            // Define color thresholds
            let blue_threshold = 0.7;
            let beach_threshold = 0.705;
            let green_threshold = 0.90;
            let brown_threshold = 0.95;

            // Define colors
            let deep_blue = Rgba([0, 0, 139, 255]);
            let light_blue = Rgba([30 as u8, 144 as u8, 255 as u8, 255]); // Light Blue

            // Assign colors based on thresholds
            let color = if height_value < blue_threshold {
                // Interpolate between deep_blue and light_blue based on height_value
                let t = (height_value - 0.3) / (0.7 - 0.3); // Adjust thresholds as needed
                interpolate_color(deep_blue, light_blue, t)
            } else if height_value < beach_threshold {
                Rgba([255, 235, 205, 255]) // Beach
            } else if height_value < green_threshold {
                Rgba([34, 139, 34, 255]) // Green
            } else if height_value < brown_threshold {
                Rgba([160, 82, 45, 255]) // Brown
            } else {
                Rgba([255, 255, 255, 255]) // White
            };

            image_buffer.put_pixel(x as u32, y as u32, color);
        }
    }

    if let Err(e) = image_buffer.save(output_path) {
        eprintln!("Error saving image: {}", e);
    }
}

pub fn generate_plane_map_builder(seed: u32, output_path: &str) {
    let fbm = Fbm::<Perlin>::new(seed); // Four octaves of FBM noise
    let ridge = RidgedMulti::<Perlin>::new(seed); // Four octaves of Ridge noise
    let fbm_file_name = format!("{}-fbm.png", output_path);
    let ridge_file_name = format!("{}-ridge.png", output_path);

    let scale = 0.5;
    PlaneMapBuilder::<&Fbm<Perlin>, 2>::new(&fbm)
        .set_size(1000, 1000)
        .set_x_bounds(-1.0 * scale, 1.0 * scale)
        .set_y_bounds(-1.0 * scale, 1.0 * scale)
        .build()
        .write_to_file(&fbm_file_name);

    PlaneMapBuilder::<&RidgedMulti<Perlin>, 2>::new(&ridge)
        .set_size(1000, 1000)
        .set_x_bounds(-1.0 * scale, 1.0 * scale)
        .set_y_bounds(-1.0 * scale, 1.0 * scale)
        .build()
        .write_to_file(&ridge_file_name);
}

pub fn write_heightmap(heightmap: &Vec<Vec<f64>>, output_path: &str) {
    let width = heightmap.len();
    let height = heightmap[0].len();

    let mut image_buffer = ImageBuffer::new(width as u32, height as u32);

    for x in 0..width {
        for y in 0..height {
            let normalized_height = heightmap[x][y];
            let color_value = (normalized_height * 255.0) as u8;
            let pixel = Rgba([color_value, color_value, color_value, 255]);

            image_buffer.put_pixel(x as u32, y as u32, pixel);
        }
    }

    if let Err(e) = image_buffer.save(output_path) {
        eprintln!("Error saving image: {}", e);
    }
}

pub fn generate_noise_map(
    map_width: usize,
    map_height: usize,
    seed: u32,
    scale: f64, // Use f64 for scale
    octaves: usize,
    persistence: f64, // Use f64 for persistence
    lacunarity: f64,
    offset: (f64, f64), // Use f64 for offsets
) -> Vec<Vec<f64>> {
    let mut noise_map = vec![vec![0.0; map_height]; map_width];

    let mut rng = rand::thread_rng();
    let mut octave_offsets = Vec::with_capacity(octaves);
    for _ in 0..octaves {
        let offset_x = rng.gen_range(-100000.0..100000.0) + offset.0;
        let offset_y = rng.gen_range(-100000.0..100000.0) + offset.1;
        octave_offsets.push((offset_x, offset_y));
    }

    if scale <= 0.0 {
        return noise_map; // Return an empty map if scale is invalid
    }

    let mut max_noise_height = f64::MIN;
    let mut min_noise_height = f64::MAX;

    let half_width = map_width as f64 / 2.0;
    let half_height = map_height as f64 / 2.0;

    let perlin = Perlin::new(seed);

    for y in 0..map_height {
        for x in 0..map_width {
            let mut amplitude = 1.0;
            let mut frequency = 1.0;
            let mut noise_height = 0.0;

            for i in 0..octaves {
                let sample_x = (x as f64 - half_width) / scale * frequency + octave_offsets[i].0;
                let sample_y = (y as f64 - half_height) / scale * frequency + octave_offsets[i].1;

                let perlin_value = perlin.get([sample_x, sample_y]) as f64; // Use f64 for the result
                noise_height += perlin_value * amplitude;

                amplitude *= persistence;
                frequency *= lacunarity;
            }

            if noise_height > max_noise_height {
                max_noise_height = noise_height;
            } else if noise_height < min_noise_height {
                min_noise_height = noise_height;
            }
            noise_map[x][y] = noise_height;
        }
    }

    for y in 0..map_height {
        for x in 0..map_width {
            noise_map[x][y] =
                (noise_map[x][y] - min_noise_height) / (max_noise_height - min_noise_height);
        }
    }

    noise_map
}
