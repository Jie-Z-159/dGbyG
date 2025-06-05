from flask import Flask
from dgbygapp.config import Config

def create_app(config_class=Config):
    app = Flask(__name__)
    app.config.from_object(config_class)
    
    from dgbygapp.routes import init_routes
    init_routes(app)
    
    return app
