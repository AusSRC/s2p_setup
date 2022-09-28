import os
import glob
import unittest
import s2p_setup


class TestSetup(unittest.TestCase):
    def setUp(self):
        """Make directory for temporary output cubelet fits files.

        """
        test_subdirs = [
            'outputs',
            'pipeline_test_outputs',
            'small_region_outputs'
        ]
        for subdir in test_subdirs:
            if os.path.exists(subdir):
                files = glob.glob(f"{subdir}/*")
                for f in files:
                    os.remove(f)
            else:
                os.mkdir(subdir)

    def tearDown(self):
        """Remove all cubelet files.

        """
        test_subdirs = [
            'outputs',
            'pipeline_test_outputs',
            'small_region_outputs'
        ]
        for subdir in test_subdirs:
            files = glob.glob(f"{subdir}/*")
            for f in files:
                os.remove(f)
            os.rmdir(subdir)

    def test_generate_config_with_region(self):
        """Run s2p_setup with region and assert that the correct sub-cube
        parameter files are created. Test that this runs to completion.
        NOTE: tested on AusSRC Carnaby with known values.

        """
        s2p_setup.main([
            "--config", '/mnt/shared/wallaby/config/s2p_setup.ini',
            "--image_cube", '/mnt/shared/wallaby/data/image.restored.i.NGC4808_A.SB33681.cube.contsub.fits',
            "--weights_cube", '/mnt/shared/wallaby/data/weights.i.NGC4808_A.SB33681.cube.fits',
            "--region", '192.4, 195.4, 4.50, 7.50',
            "--run_name", 'TestSetup',
            "--sofia_template", '/mnt/shared/wallaby/config/sofia.par',
            "--output_dir", 'outputs',
            "--products_dir", 'products',
        ])

    def test_number_of_subcubes(self):
        """Two part test.
        1. Test behaviour where no --region provided (default in config is the entire cube)
        2. Assert number of parameter files generated for the region is as expected.

        N_expect determined for sub-cube size in configuration 1500 1500 1400

        """
        N_expect = 54
        s2p_setup.main([
            "--config", '/mnt/shared/wallaby/config/s2p_setup.ini',
            "--image_cube", '/mnt/shared/wallaby/data/image.restored.i.NGC4808_A.SB33681.cube.contsub.fits',
            "--weights_cube", '/mnt/shared/wallaby/data/weights.i.NGC4808_A.SB33681.cube.fits',
            "--run_name", 'TestSetup',
            "--sofia_template", '/mnt/shared/wallaby/config/sofia.par',
            "--output_dir", 'outputs',
            "--products_dir", 'products',
        ])
        parameter_files = glob.glob("outputs/sofia_*.par")
        self.assertEqual(len(parameter_files), N_expect)

    def test_region_smaller_than_subcube(self):
        """Test behaviour where the region is smaller in each dimension than the minimum required sub-cube size.
        Expect 1 x 1 x 6 regions.

        """
        N_exp = 6
        s2p_setup.main([
            "--config", '/mnt/shared/wallaby/config/s2p_setup.ini',
            "--image_cube", '/mnt/shared/wallaby/data/image.restored.i.NGC4808_A.SB33681.cube.contsub.fits',
            "--weights_cube", '/mnt/shared/wallaby/data/weights.i.NGC4808_A.SB33681.cube.fits',
            "--region", '195.0, 196.0, 5.0, 6.0',
            "--run_name", 'TestSetup',
            "--sofia_template", '/mnt/shared/wallaby/config/sofia.par',
            "--output_dir", 'small_region_outputs',
            "--products_dir", 'products',
        ])
        parameter_files = glob.glob("small_region_outputs/sofia_*.par")
        self.assertEqual(len(parameter_files), N_exp)

    def test_WALLABY_pipeline_test_case(self):
        """Assert this code works for the WALLABY pipeline test-case in the AusSRC system.

        """
        s2p_setup.main([
            "--config", '/mnt/shared/wallaby/config/s2p_setup.milkyway.ini',
            "--region", '200.77, 204.77, -24.605, -20.605',
            "--image_cube", '/mnt/shared/wallaby/post-runs/test/mosaic.fits',
            "--weights_cube", '/mnt/shared/wallaby/post-runs/test/weights.mosaic.fits',
            "--run_name", 'TestSetup',
            "--sofia_template", '/mnt/shared/wallaby/config/sofia.par',
            "--output_dir", 'pipeline_test_outputs',
            "--products_dir", 'products',
        ])
