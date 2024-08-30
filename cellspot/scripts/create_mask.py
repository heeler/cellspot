from cellspot.spot_finder import SpotFinder

from pathlib import Path
import argparse


def main(namespace_args):


    sf = SpotFinder(
        image_path=namespace_args.filename
    )
    mask = sf.run(clip_values=namespace_args.clip_to,
                  red_threshold=namespace_args.red_threshold,
                  pixel_cycles=namespace_args.pixel_cycles)

    output = (namespace_args.output_folder / namespace_args.filename.name).with_suffix('.json')
    return sf.write_results(file_out=output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='find nuclei with associated inclusion bodies')
    parser.add_argument('filename', action='store', type=Path)
    parser.add_argument('-r', '--red_threshold', action='store', default=82, type=int)
    parser.add_argument('-c', '--clip_to', action='store', default=(5.0, 90.0), type=tuple[float, float])
    parser.add_argument('-o', '--output_folder', action='store', default=False, type=Path, required=True)
    parser.add_argument('-p', '--pixel_cycles', action='store', default=6, type=int)
    args = parser.parse_args()
    main(args)
