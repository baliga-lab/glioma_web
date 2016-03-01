import sys, os
sys.path.insert(0, '/local/apache-stuff/glioma_web')
os.environ['PYTHON_EGG_CACHE'] = '/local/apache-stuff/python_egg_cache'
os.environ['GLIOMA_SETTINGS'] = '/local/apache-stuff/glioma_web_stati/settings.cfg'
from app import app as application
application.debug = True
application.secret_key = 'trsntsrotrsenoiwfpguy8fpgfp'
