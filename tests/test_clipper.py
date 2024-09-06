#!/usr/bin/env python

# import os
import unittest
import bobkit


class TestClipper(unittest.TestCase):
    def test_clipper_core(self):
        c = bobkit.clipper.Test_core()
        self.assertTrue(c())

    def test_clipper_contrib(self):
        c = bobkit.clipper.Test_contrib()
        self.assertTrue(c())

    # def test_clipper_minimol_gemmi(self):
    #    c = bobkit.clipper.Test_minimol_gemmi()
    #    result = c.run("")
    #    self.assertTrue()


if __name__ == "__main__":
    unittest.main()
