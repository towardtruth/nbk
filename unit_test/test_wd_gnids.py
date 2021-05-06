import unittest
from Data import GetData as GD
import settings


class TestWDGnids(unittest.TestCase):
    def test_wd_all_gnids(self):
        settings.conn_string = settings.conn_string_test
        all_gnids = GD.wd_all_gnid(gene_prot='g')
        print(len(all_gnids))


if __name__ == '__main__':
    unittest.main()
