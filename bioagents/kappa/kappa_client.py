"""Web API client for a Kappa simulator."""

import urllib, urllib2
import json
import requests
from time import sleep

class RuntimeError(Exception):
    def __init__(self, errors):
        self.errors = errors

kappa_default = 'http://api.executableknowledge.org/kappa/v2/projects/default'

class KappaRuntime(object):
    def __init__(self, endpoint=None):
        """Create a Kappa client."""
        if not endpoint:
            self.url = kappa_default
        else:
            self.url = endpoint

    def version(self):
        """Return the version of the Kappa environment."""
        try:
            res = requests.get(self.url)
            if res.status_code != 200:
                raise Exception('Kappa service returned with code: %s' %
                                res.status_code)
            content = res.json()
            version = content.get('project_version')
            return version
        except Exception as e:
            raise RuntimeError(e)

    def parse(self, code):
        """Parse given Kappa model code and throw exception if fails."""
        query_args = {'code': code}
        parse_url = self.url + '/parse'
        res = requests.post(parse_url, json=query_args)
        #try:
        content = response.json()
        return json.loads(content)
        #except urllib2.HTTPError as e:
        #    if e.code == 400:
        #        error_details = json.loads(e.read())
        #        raise RuntimeError(error_details)
        #    else:
        #        raise e
        #except urllib2.URLError as e:
        #    RuntimeError(e.reason)

    def start(self, parameter):
        """Start a simulation with given parameters."""
        if not 'max_time' in parameter:
            parameter['max_time'] = None
        if not 'max_events' in parameter:
            parameter['max_events'] = None
        code = json.dumps(parameter)
        method = "POST"
        handler = urllib2.HTTPHandler()
        opener = urllib2.build_opener(handler)
        parse_url = "{0}/process".format(self.url)
        request = urllib2.Request(parse_url, data=code)
        request.get_method = lambda: method
        try:
            connection = opener.open(request)
        except urllib2.HTTPError,e:
            connection = e
        except urllib2.URLError as e:
            raise RuntimeError(e.reason)

        if connection.code == 200:
            text = connection.read()
            return int(json.loads(text))
        elif connection.code == 400:
            text = connection.read()
            error_details = json.loads(text)
            raise RuntimeError(error_details)
        else:
            raise e

    def stop(self, token):
        """Stop a given simulation."""
        method = "DELETE"
        handler = urllib2.HTTPHandler()
        opener = urllib2.build_opener(handler)
        parse_url = "{0}/process/{1}".format(self.url,token)
        request = urllib2.Request(parse_url)
        request.get_method = lambda: method
        try:
            connection = opener.open(request)
        except urllib2.HTTPError,e:
            connection = e
        except urllib2.URLError as e:
            raise RuntimeError(e.reason)

        if connection.code == 200:
            text = connection.read()
            return None
        elif connection.code == 400:
            text = connection.read()
            error_details = json.loads(text)
            raise RuntimeError(error_details)
        else:
            raise e

    def status(self, token):
        """Return status of running simulation."""
        try:
            version_url = "{0}/process/{1}".format(self.url,token)
            response = urllib2.urlopen(version_url)
            text = response.read()
            return json.loads(text)
        except urllib2.HTTPError as e:
            if e.code == 400:
                error_details = json.loads(e.read())
                raise RuntimeError(error_details)
            else:
                raise e
        except urllib2.URLError as e:
            RuntimeError(e.reason)


    def shutdown(self,key):
        """Shutdown the server."""
        method = "POST"
        handler = urllib2.HTTPHandler()
        opener = urllib2.build_opener(handler)
        parse_url = "{0}/shutdown".format(self.url)
        request = urllib2.Request(parse_url, data=key)
        request.get_method = lambda: method
        try:
            connection = opener.open(request)
        except urllib2.HTTPError,e:
            connection = e
        except urllib2.URLError as e:
            raise RuntimeError(e.reason)
        if connection.code == 200:
            text = connection.read()
            return text
        elif connection.code == 400:
            text = connection.read()
            raise RuntimeError(text)
        elif connection.code == 401:
            text = connection.read()
            raise RuntimeError(text)
        else:
            raise e

if __name__ == "__main__":
    with open("../abc-pert.ka") as f:
        try:
            subprocess.Popen('../WebSim.native --shutdown-key 6666 --port 6666'.split())
            sleep(1)
            data = f.read()
            runtime = KappaRuntime("http://localhost:6666")
            token = runtime.start({ 'code': data
                                  , 'nb_plot': 10
                                  , 'max_events' : 10000 })
            sleep(10)
            status = runtime.status(token)
            print status
            #print render_status(status).toString()
            print runtime.shutdown('6666')
        except RuntimeError as e:
            print e.errors
