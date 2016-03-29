import sys, os
os.environ['R_HOME'] = '/tools/R-3.0.3/lib64/R'
sys.path.insert(0, '/local/apache-stuff/glioma_web')
os.environ['PYTHON_EGG_CACHE'] = '/local/apache-stuff/python_egg_cache'
os.environ['GLIOMA_SETTINGS'] = '/local/apache-stuff/glioma_web_data/settings.cfg'
from app import app as application
application.debug = True
application.secret_key = 'trsntsrotrsenoiwfpguy8fstrtpgfp'
