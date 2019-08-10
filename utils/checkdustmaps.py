"""
If the SFD or Bayestar maps haven't been downloaded, then this script downloads them.
"""
import os

__all__ = ['check_dustmaps']

def check_dustmaps(_path):
    from dustmaps.config import config
    print(_path)
    config['data_dir'] = _path

    if not os.path.exists(config['data_dir']+'sfd'):
        import dustmaps.sfd
        dustmaps.sfd.fetch()

    #if not os.path.exists(config['data_dir']+'bayestar'):
    #    import dustmaps.bayestar
    #    dustmaps.bayestar.fetch()

