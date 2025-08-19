#!/usr/bin/env python3
"""
Tissuumaps Manager Utility
Handles Tissuumaps server management and integration
"""

import streamlit as st
import subprocess
import time
import webbrowser
from pathlib import Path

class TissuumapsManager:
    """Manages Tissuumaps integration"""
    
    def __init__(self, data_dir: str):
        self.data_dir = Path(data_dir)
        self.server_process = None
        self.port = 5100
        
    def start_tissuumaps_server(self, port: int = None) -> bool:
        """Start Tissuumaps server"""
        try:
            if port:
                self.port = port
            
            # Kill any existing process on the port
            subprocess.run(f"lsof -ti:{self.port} | xargs kill -9", shell=True, capture_output=True)
            
            # Start Tissuumaps server
            cmd = f"tissuumaps --port {self.port} --path {self.data_dir}"
            self.server_process = subprocess.Popen(cmd, shell=True)
            
            # Wait a moment for server to start
            time.sleep(3)
            
            return True
        except Exception as e:
            st.error(f"Error starting Tissuumaps server: {e}")
            return False
    
    def stop_server(self):
        """Stop Tissuumaps server"""
        if self.server_process:
            self.server_process.terminate()
            self.server_process = None
    
    def get_tissuumaps_url(self, port: int = None) -> str:
        """Get Tissuumaps URL"""
        if port:
            self.port = port
        return f"http://localhost:{self.port}"
    
    def is_server_running(self) -> bool:
        """Check if Tissuumaps server is running"""
        try:
            # Check if process is still running
            if self.server_process:
                return self.server_process.poll() is None
            return False
        except Exception:
            return False
    
    def get_server_status(self) -> dict:
        """Get server status information"""
        status = {
            'running': self.is_server_running(),
            'port': self.port,
            'url': self.get_tissuumaps_url(),
            'data_dir': str(self.data_dir)
        }
        return status
