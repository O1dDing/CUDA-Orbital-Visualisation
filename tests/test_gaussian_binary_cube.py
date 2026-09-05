from pathlib import Path
import struct
import tempfile
import unittest

from gaussian_binary_cube import np, read_binary_mo_cube, write_single_orbital_cube
from grid_reference_compare import read_mo_cube


def packet(payload):
    marker = struct.pack('<I', len(payload))
    return marker + payload + marker


def specimen(*, duplicate_ids=False, nonfinite=False, nval=2, values=None):
    chunks = [packet(text.encode().ljust(80)) for text in ('Known binary fixture', 'Alpha MO coefficients')]
    chunks.append(packet(struct.pack('<q3dq', -1, -1., -2., -3., nval)))
    for row in ((2, .5, 0., 0.), (2, 0., .5, 0.), (2, 0., 0., .5)):
        chunks.append(packet(struct.pack('<q3d', *row)))
    chunks.append(packet(struct.pack('<q4d', 1, 1., 0., 0., 0.)))
    chunks.append(packet(struct.pack('<3q', 2, 7, 7 if duplicate_ids else 12)))
    if values is None:
        values = np.column_stack((np.arange(1., 9.), np.arange(11., 19.)))
    if nonfinite:
        values[3, 1] = np.nan
    for row in range(4):
        chunks.append(packet(values[row*2:row*2+2].astype('<f8').tobytes()))
    return b''.join(chunks)


class GaussianBinaryCubeTests(unittest.TestCase):
    def read(self, payload):
        with tempfile.TemporaryDirectory() as raw:
            path = Path(raw) / 'cube.cub'
            path.write_bytes(payload)
            return read_binary_mo_cube(path)

    def test_binary_rows_preserve_point_dataset_interleaving(self):
        cube = self.read(specimen())
        self.assertEqual(cube['orbital_ids'], [7, 12])
        self.assertEqual(cube['shape'], [2, 2, 2])
        np.testing.assert_array_equal(cube['values'], [np.arange(1., 9.), np.arange(11., 19.)])

    def test_malformed_markers_counts_and_values_cannot_be_accepted(self):
        complete = specimen()
        for payload in (complete[:-1], complete+b'X', b'\x00\x00\x00\x50'+complete[4:],
                        complete[:84]+b'\x01\x00\x00\x00'+complete[88:],
                        specimen(duplicate_ids=True), specimen(nonfinite=True), specimen(nval=3)):
            with self.subTest(length=len(payload)), self.assertRaises(ValueError):
                self.read(payload)

    def test_memory_mapped_grid_is_exact_and_never_overwrites(self):
        with tempfile.TemporaryDirectory() as raw:
            folder = Path(raw)
            cube = folder / 'cube.cub'
            cube.write_bytes(specimen())
            target = folder / 'values.npy'
            result = read_binary_mo_cube(cube, values_path=target)
            np.testing.assert_array_equal(np.load(target), result['values'])
            with self.assertRaises(FileExistsError):
                read_binary_mo_cube(cube, values_path=target)
            result['values']._mmap.close()

    def test_single_dataset_packaging_retains_exact_float64_bits(self):
        with tempfile.TemporaryDirectory() as raw:
            folder = Path(raw)
            source = folder / 'all.cub'
            values = np.column_stack((np.arange(1., 9.),
                [np.nextafter(1., 2.), -0., 1e-250, -1e250, np.pi, -np.pi, 1./3, 7.]))
            source.write_bytes(specimen(values=values))
            cube = read_binary_mo_cube(source)
            target = write_single_orbital_cube(cube, 1, folder / 'one.cube')
            text = read_mo_cube(target)
            self.assertEqual(text['orbital_ids'], [12])
            np.testing.assert_array_equal(text['values'][0].view(np.uint64), cube['values'][1].view(np.uint64))
            with self.assertRaises(FileExistsError):
                write_single_orbital_cube(cube, 1, target)


if __name__ == '__main__':
    unittest.main()
