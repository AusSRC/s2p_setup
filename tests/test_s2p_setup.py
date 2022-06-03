import os
import glob
import unittest
import s2p_setup


class TestSetup(unittest.TestCase):
    def setUp(self):
        """Make directory for temporary output cubelet fits files.

        """
        if os.path.exists("outputs"):
            files = glob.glob("outputs/*")
            for f in files:
                os.remove(f)
        else:
            os.mkdir("outputs")

    def tearDown(self):
        """Remove all cubelet files.

        """
        files = glob.glob("outputs/*")
        for f in files:
            os.remove(f)
        os.rmdir("outputs")

    def test_generate_config_with_region(self):
        """Run s2p_setup with region and assert that the correct sub-cube
        parameter files are created. Test that this runs to completion.
        NOTE: tested on AusSRC Carnaby with known values.

        """
        s2p_setup.main([
            "--config", '/mnt/shared/wallaby/config/s2p_setup.ini',
            "--image_cube", '/mnt/shared/wallaby/data/image.restored.i.NGC4808_A.SB33681.cube.contsub.fits',
            "--region", '194.4, 195.4, 5.50, 6.50',
            "--run_name", 'TestSetup',
            "--sofia_template", '/mnt/shared/wallaby/config/sofia.par',
            "--output_dir", 'outputs',
            "--products_dir", 'products',
        ])

    # TODO(austin): Test functionality where region not provided (should just use pre-set pixel coordinates)