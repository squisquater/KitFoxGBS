import rasterio
import matplotlib.pyplot as plt
import numpy as np
import argparse

def plot_raster(input_tif, output_png):
    with rasterio.open(input_tif) as src:
        data = src.read(1)

        # Check the data range
        print(f"Data min: {data.min()}, Data max: {data.max()}")

        # Plotting
        plt.figure(figsize=(10, 10))
        plt.imshow(data, cmap='gray', vmin=data.min(), vmax=data.max())
        plt.colorbar()
        plt.title('California Major Roads Rasterized')
        plt.xlabel('Easting')
        plt.ylabel('Northing')

        # Save the plot
        plt.savefig(output_png)
        plt.show()

def main():
    parser = argparse.ArgumentParser(description="Plot and save rasterized TIFF file as PNG")
    parser.add_argument('input_tif', type=str, help='Path to the input .tif file')
    parser.add_argument('output_png', type=str, help='Path to the output .png file')

    args = parser.parse_args()
    plot_raster(args.input_tif, args.output_png)

if __name__ == '__main__':
    main()
