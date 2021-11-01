#!/usr/bin/env python
# $Id: iprscan5_soappy.py 2760 2014-04-10 15:24:31Z hpm $
# ======================================================================
#
# Copyright 2008-2014 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# ======================================================================
# InterProScan 5 (SOAP) service, Python client using SOAPpy.
#
# Tested with:
#   Python 2.5.2 with SOAPpy 0.12.0 (Ubuntu 8.04 LTS)
#   Python 2.6.5 with SOAPpy 0.12.0 (Ubuntu 10.04 LTS)
#   Python 2.7.3 with SOAPpy 0.12.0 (Ubuntu 12.04 LTS)
#
# See:
# http://www.ebi.ac.uk/Tools/webservices/services/pfa/iprscan5_soap
# http://www.ebi.ac.uk/Tools/webservices/tutorials/06_programming/python
# ======================================================================
# WSDL URL for service
# edited by niek de klein
wsdlUrl = 'http://www.ebi.ac.uk/Tools/services/soap/iprscan5?wsdl'

# Load libraries
import base64, platform, os, SOAPpy, time
from SOAPpy import WSDL
import urllib2
import xml.parsers.expat
import sys
class InterproScan():
    """Identify protein family, domain and signal signatures in a 
       protein sequence using InterProScan. For more information on InterProScan 
       refer to http://www.ebi.ac.uk/Tools/pfa/iprscan
       
       For further information about the InterProScan (SOAP) web service, see http://www.ebi.ac.uk/Tools/webservices/services/pfa/iprscan_soap.
    
       Args:
           checkInterval (int)       Set interval for checking status (Default = 3)
           outputLevel (int)         Output level (Default = 1)
           self.debugLevel (int)          Debug level (Default = 0)
       
           # Tool specific options
           appl (list)               signature methods to use, see paramDetail = 'appl'
           crc  (bool)               enable InterProScan Matches look-up (faster)
           nocrc (bool)              disable InterProScan Matches look-up (slower)
           goterms (bool)            enable inclusion of GO terms
           nogoterms (bool)          disable inclusion of GO terms
           sequence_file (str)       input sequence file name
           # General options
           email (str)               e-mail address
           title (str)               job title
           outfile (str)             file name for results
           outformat (str)           output format for results
           async (bool)              asynchronous mode
           jobid (str)               job identifier
           polljob (bool)            get job result
           status (bool)             get job status
           resultTypes (bool)        get result types
           params (bool)             list input parameters
           paramDetail (str)         get details for parameter
           quiet (bool)              decrease output level
           verbose (bool)            increase output level
           trace (bool)              show SOAP messages
           sequences (str)       List of sequences to search interpro with. 
           WSDLUrl (str)             WSDL URL for service (Default = 'http://www.ebi.ac.uk/Tools/services/soap/iprscan5?wsdl')
           self.debugLevel (int)          debug output level
           write_outfiles (bool)     Write the result files out

    """
    
    def __init__(self, checkInterval = 3, outputLevel = 1, debugLevel = 1, appl = None, crc = False, nocrc = False,
                 goterms = False, nogoterms = False, sequence_file = None, email = None, title = None, outfile = None, 
                 outformat = None, async = True, jobid = None, polljob = False, status = False, resultTypes = False, 
                 params = False, paramDetail = None, quiet = False, verbose = False, trace = False, sequences = None,
                 WSDLUrl = 'http://www.ebi.ac.uk/Tools/services/soap/iprscan5?wsdl', write_outfiles = True):
        self.outfile = outfile
        self.outformat = outformat
        self.debugLevel = debugLevel
        self.checkInterval = checkInterval
        self.outputLevel = outputLevel
        self.sequence_name_for_job_id = {}
        # Increase output level
        if verbose:
            self.outputLevel += 1
        # Decrease output level
        if quiet:
            self.outputLevel -= 1
        # Set the client user-agent.
        clientRevision = '$Revision: 2760 $'
        clientVersion = '0'
        if len(clientRevision) > 11:
            clientVersion = clientRevision[11:-2]
        self.userAgent = 'EBI-Sample-Client/%s (%s; Python %s; %s)' % (
                                                             clientVersion, os.path.basename( __file__ ),
                                                             platform.python_version(), platform.system()
                                                             )
        # Redefine default User-agent function to return custom User-agent.
        SOAPpy.Client.SOAPUserAgent = self.SOAPUserAgent
        self.printDebugMessage('main', 'User-agent: ' + SOAPpy.Client.SOAPUserAgent(), 1)
        
        # Create the service interface
        self.printDebugMessage('main', 'WSDL: ' + WSDLUrl, 1)
        getServer = True
        while getServer:
            try:
                self.server = WSDL.Proxy(WSDLUrl)
                getServer = False
            except xml.parsers.expat.ExpatError:
                print('WSDL not well formed, possibly no internet connection. Trying again in 5 seconds')
                time.sleep(5)
        # Fix message namespace (not set from the WSDL).
        for method in self.server.methods:
            if self.server.methods[method].namespace == None:
                self.server.methods[method].namespace = 'http://soap.jdispatcher.ebi.ac.uk'
        
        # Configure HTTP proxy from OS environment (e.g. http_proxy="http://proxy.example.com:8080")
        if os.environ.has_key('http_proxy'):
            http_proxy_conf = os.environ['http_proxy'].replace('http://', '')
        elif os.environ.has_key('HTTP_PROXY'):
            http_proxy_conf = os.environ['HTTP_PROXY'].replace('http://', '')
        else:
            http_proxy_conf = None
        self.server.soapproxy.http_proxy = http_proxy_conf
        
        # If required enable SOAP message trace
        if trace:
            self.server.soapproxy.config.dumpSOAPOut = 1
            self.server.soapproxy.config.dumpSOAPIn = 1
        
        # List parameters
        if params:
            for paramName in self.serviceGetParameters()['id']:
                sys.stderr.write(str(paramName)+'\n')
        # Get parameter details
        elif paramDetail:
            self.printGetParameterDetails(paramDetail)
            exit(-1)
        # Submit job
        elif email and not jobid:
            params = {}
            if not sequence_file and not sequences:
                raise ValueError('sequence_file and single_sequence where empty. One of them needs input')
            # get list of sequences from the file
            if sequence_file:
                sequences = self.readFile(sequence_file)
            elif isinstance(sequences, list) == False:
                sequences = [sequences]
            # Booleans need to be represented as 1/0 rather than True/False
            if crc:
                params['nocrc'] = 0
            if nocrc:
                params['nocrc'] = 1
            if goterms:
                params['goterms'] = 1
            if nogoterms:
                params['goterms'] = 0
            # Add the other options (if defined)
            if appl:
                params['appl'] = {'string':appl}
            # Submit the job
            self.jobid_list = []
            self.sequence_for_job = {}
            for sequence in sequences:
                params['sequence'] = sequence
                jobid = self.serviceRun(email, title, params)
                self.sequence_name_for_job_id[jobid] = sequence.split('\n')[0].lstrip('>')
                self.sequence_for_job[jobid] = sequence
                self.jobid = jobid
                if self.outputLevel > 0:
                    sys.stderr.write(str(jobid)+'\n')
                self.jobid_list.append(jobid)
                if write_outfiles:
                    self.getResult(jobid)

            # check if all the jobs executed correctly, otherwise run again
            jobid_new = None
            for jobid in self.jobid_list:
                error = True
                while error:
                    try:
                        if jobid_new != None:
                            self.sequence_for_job[jobid_new] = self.sequence_for_job[jobid]
                            self.jobid_list.append(jobid_new)
                            self.jobid_list.remove(jobid)
                            jobid = jobid_new
                        if not self.clientPoll(jobid) == 'ERROR':    # wait till job is finished     
                            # check if xml is in the type format, otherwise have to run it again
                            for result_type in self.serviceGetResultTypes(jobid):
                                try:
                                    if result_type['identifier'] == 'xml':
                                        error = False
                                        break   # if it works, break out of the while 1 loop
                                except TypeError:
                                    sys.stderr.write('Error with result_type, try again later. Error happened with sequence:\n'+str(sequence)+'\n')
                            if error == False:
                                break
                            sys.stderr.write('xml was not found in serviceGetResultTypes\n')
                    except xml.parsers.expat.ExpatError:                                     # if this error is thrown
                        sys.stderr.write('parser error\n')                  #    print that the error was thrown
                    sys.stderr.write('Trying again after 1 seconds\n')
                    time.sleep(1)   # per time it goes wrong, wait 10 sec longer, so it doesn't keep spamming loops
                    if self.sequence_for_job.has_key(jobid):
                        params['sequence'] = self.sequence_for_job[jobid]
                    else:
                        jobid_new = jobid
                    error = True
                    jobid_new = self.serviceRun(email, title, params)
                jobid_new = None
        # Get job status
        elif status and jobid:
            sys.stderr.write('Checking status')
            status = self.serviceCheckStatus(jobid)
            sys.stderr.write(str(status)+'\n')
            time.sleep(5)
        # List result types for job
        elif resultTypes and jobid:
            self.printGetResultTypes(jobid)
        # Get results for job
        elif polljob and jobid:
            if write_outfiles:
                self.getResult(jobid)
        else:
            sys.stderr.write('Error: unrecognised argument combination\n')
    
    def getSequenceNameForJobID(self,jobid):
        if jobid in self.sequence_name_for_job_id:
            return(self.sequence_name_for_job_id[jobid])
        else:
            print('jobid '+jobid+' not found')
    
    def getFastaSequenceForJobID(self, jobid):
        if jobid in self.sequence_for_job:
            return(self.sequence_for_job[jobid])
        else:
            return None
    
    def getJobIDs(self):
        return self.jobid_list
    
    def SOAPUserAgent(self):
        return self.userAgent
    # Debug print
    def printDebugMessage(self,functionName, message, level):
        if(level <= self.debugLevel):
            print >>sys.stderr, '[' + functionName + '] ' + message

    # Get input parameters list
    def serviceGetParameters(self):
        self.printDebugMessage('serviceGetParameters', 'Begin', 1)
        result = self.server.getParameters()
        self.printDebugMessage('serviceGetParameters', 'End', 1)
        return result

    # Get input parameter information
    def serviceGetParameterDetails(self,paramName):
        self.printDebugMessage('serviceGetParameterDetails', 'Begin', 1)
        result= self.server.getParameterDetails(parameterId=paramName)
        self.printDebugMessage('serviceGetParameterDetails', 'End', 1)
        return result

    # Submit job
    def serviceRun(self, email, title, params):
        self.printDebugMessage('serviceRun', 'Begin', 1)
        run = True
        while run:
            try:
                jobid = self.server.run(email=email, title=title, parameters=params)
                run = False
            except:
                print('job cannot be run on server. Possibly no internet connection or Interpro server is not responding. Trying again in 5 seconds')
                time.sleep(5)
    
        self.printDebugMessage('serviceRun', 'End', 1)
        return jobid

    # Get job status
    def serviceCheckStatus(self, jobId):
        self.printDebugMessage('serviceCheckStatus', 'jobId: ' + jobId, 1)
        try:
            try:
                result = self.server.getStatus(jobId = jobId)
            except:
                print('sleeping 30 seconds before checking status again...')
                time.sleep(30)
                result = self.server.getStatus(jobId = jobId)
        except:
            result = None
        return result

    # Get available result types for job
    def serviceGetResultTypes(self, jobId):
        self.printDebugMessage('serviceGetResultTypes', 'Begin', 1)
        result = self.server.getResultTypes(jobId=jobId)
        self.printDebugMessage('serviceGetResultTypes', 'End', 1)
        return result['type']

    # Get result
    def serviceGetResult(self, jobId, service_type):
        self.printDebugMessage('serviceGetResult', 'Begin', 1)
        self.printDebugMessage('serviceGetResult', 'jobId: ' + jobId, 1)
        self.printDebugMessage('serviceGetResult', 'type: ' + service_type, 1)
        self.resultBase64 = self.server.getResult(jobId=jobId, type=service_type)
        result = base64.decodestring(self.resultBase64)
        self.printDebugMessage('serviceGetResult', 'End', 1)
        return result

    # Client-side poll
    def clientPoll(self, jobId):
        self.printDebugMessage('clientPoll', 'Begin', 1)
        result = 'PENDING'
        while result == 'RUNNING' or result == 'PENDING':
            sys.stderr.write('Checking status\n')
            result = self.serviceCheckStatus(jobId)
            print >>sys.stderr, result
            if result == 'RUNNING' or result == 'PENDING':
                print('sleeping 30 seconds before checking status again...')
                time.sleep(30)
        self.printDebugMessage('clientPoll', 'End', 1)

    # Get result for a jobid
    def getResult(self, jobId):
        self.printDebugMessage('getResult', 'Begin', 1)
        self.printDebugMessage('getResult', 'jobId: ' + jobId, 1)
        # Check status and wait if necessary
        self.clientPoll(jobId)
        # Get available result types
        resultTypes = self.serviceGetResultTypes(jobId)
        for resultType in resultTypes:
            # Get the result
            result = self.serviceGetResult(jobId, resultType['identifier'])
            # Derive the filename for the result
            # Write a result file
        self.printDebugMessage('getResult', 'End', 1)
        return result
    # Read a file
    def readFile(self, filename):
        self.printDebugMessage('readFile', 'Begin', 1)
        fh = open(filename, 'r')
        data = fh.read()
        fh.close()
        self.printDebugMessage('readFile', 'End', 1)
        return data

    # Output parameter details.
    def printGetParameterDetails(self, paramName):
        self.printDebugMessage('printGetParameterDetails', 'Begin', 1)
        paramDetail = self.serviceGetParameterDetails(paramName)
        print paramDetail['name'], "\t", paramDetail['type']
        print paramDetail['description']
        for value in paramDetail['values']['value']:
            print value['value'],
            if(value['defaultValue'] == 'true'):
                print '(default)',
            print
            print "\t", value['label']
            if(hasattr(value, 'properties')):
                if(isinstance(value['properties']['property'], (list, tuple))):
                    for wsProperty in value['properties']['property']:
                        print "\t", wsProperty['key'], "\t", wsProperty['value']
                else:
                    print "\t", value['properties']['property']['key'], "\t", value['properties']['property']['value']
        self.printDebugMessage('printGetParameterDetails', 'End', 1)

    # Output available result types for job.
    def printGetResultTypes(self, jobId):
        self.printDebugMessage('printGetResultTypes', 'Begin', 1)
        for resultType in self.serviceGetResultTypes(jobId):
            print resultType['identifier']
            if(hasattr(resultType, 'label')):
                print "\t", resultType['label']
            if(hasattr(resultType, 'description')):
                print "\t", resultType['description']
            if(hasattr(resultType, 'mediaType')):
                print "\t", resultType['mediaType']
            if(hasattr(resultType, 'fileSuffix')):
                print "\t", resultType['fileSuffix']
        self.printDebugMessage('printGetResultTypes', 'End', 1)

def runInterpro(sequence_list, email):
    '''Runs one sequence through interproscan, and returns the result
    
       Args:
           sequences (str):    Fasta representation of a list of sequences
           email (string):            email adres that will be send to interproscan
           
       Returns:
           Dictionary with as key sequence name and as value a dictionary with IPR id and value a dictionary with as key dbname and as value protein name
    '''    
    internet_error = True
    while internet_error:                               # need to keep trying InterproScan until it finishes without error, because of a rare error that can happen for a reason I don't know (xml from InterPro scan not well formed)
        # check if there is an internet connecion
        try:
            urllib2.urlopen('http://google.com',timeout=1) # http://74.125.113.99 is googles ip
            internet_error = False
        except (urllib2.URLError, urllib2.HTTPError): 
            sys.stderr.write('No internet\n')
            # wait 10 seconds * loop_count, so that if the error is thrown quickly it doesn't try too often then try again
            sys.stderr.write('Trying again after:10 seconds\n')
            time.sleep(10)   # per time it goes wrong, wait 10 sec longer, so it doesn't keep spamming loops
    interpro_scan_error = True
    while interpro_scan_error:
        try:
            interproscan = InterproScan(sequences = sequence_list, email=email,  # run interpro scan with one sequence (has to be fasta format). Don't
                                                               verbose = True, write_outfiles = False, quiet = False, debugLevel = 1)     # write out any of the results, because it gets parsed directly
            interpro_scan_error = False
        except AttributeError as e:
            sys.stderr.write(str(e)+'\n')
            sys.stderr.write('Trying again after 1 second\n')
            time.sleep(1)
    for jobid in interproscan.getJobIDs():
        error = True
        while error:
            try:
                result = interproscan.serviceGetResult(jobid, 'tsv')       # get the xml result
                error = False
            except Exception as e:
                sys.stderr.write(str(e)+'\n')
                sys.stderr.write('waiting 1 sec, trying again\n')
                time.sleep(1)
                seq = interproscan.getFastaSequenceForJobID(jobid)
                if not seq:
                    result = None
                    break
                interproscan = InterproScan(sequences = interproscan.getFastaSequenceForJobID(jobid), email='spam@gmail.com',  # run interpro scan with one sequence (has to be fasta format). Don't
                                                               verbose = True, write_outfiles = False, quiet = False, debugLevel = 1)
        if not result:
            continue
        interpro_result = {}
        for domain_info in result.split('\n'):
            domain_info = domain_info.split('\t')
            if len(domain_info) > 4:
                try:
                    ipr = domain_info[11].lower()# get the ipr id
                except IndexError:
                    ipr = 'noIPR'
                ipr_name = domain_info[5].lower()
                if ipr != 'noIPR':
                    if ipr_name != '':
                        if ipr in interpro_result and 'ipr_names' in interpro_result[ipr]:
                                interpro_result[ipr]['ipr_names'].append(ipr_name.lower())
                        else:
                            interpro_result[ipr] = {'ipr_names':[ipr_name.lower()]}
                    dbname = domain_info[3]                       # get the database name
                    domain_name = domain_info[12].lower()
                    if domain_name != '':
                        if ipr in interpro_result:
                            if dbname in interpro_result[ipr]:
                                interpro_result[ipr][dbname].append(domain_name.lower())   # update with databasename and domain name
                            else:
                                interpro_result[ipr].update({dbname:[domain_name.lower()]})   # update with databasename and domain name
                        else:
                            interpro_result[ipr] = {dbname:[domain_name.lower()]}   # update with databasename and domain name

        yield interproscan.getSequenceNameForJobID(jobid), interpro_result
