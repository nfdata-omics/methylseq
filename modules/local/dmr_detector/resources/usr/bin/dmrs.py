#!/usr/bin/env python3

import sys
import csv
import statistics as stats
import argparse

class DMRfinder():
    # These attributes are common to all DMR finders
    outfile       = ""          # Name of output file
    outstream     = None        # Output stream
    current_sites = []          # Sites in current DMR
    current_chrom = ""          # Chromosome being looked at
    
    def open(self):
        self.outstream = open(self.outfile, "w")

    def close(self):
        self.outstream.close()

    def add_site(self, site):
        """
        Site is a list: [chrom, start, end, dm]
        """
        pass

    def report_dmr(self):
        """
        Write current DMR to output file.
        """
        pass

    def reset(self):
        """
        Clear DMR information to start over.
        """
        pass
    
class StreamingDMRfinder(DMRfinder):
    in_dmr    = False     # True if we are inside a putative DMR
    direction = 1         # 1 for hypermethylation, -1 for hypomethylation
    threshold = 0.3       # sites with abs(diffmeth) > threshold are considered good
    nopen     = 3         # number of sites required to open a DMR
    nskip     = 1         # close DMR if more than nskip consecutive bad sites are found
    maxdist   = 0         # if non-zero, maximum distance between two consecutive good sites
    symmetrical = False   # if True, require nopen good sites at the end of DMR too
    slop      = 0         # add bases before and after DMR (fixed number if > 1, as fraction if < 1)
    
    _ngood    = 0
    _nbad     = 0

    def reset(self):
        self.current_sites = []
        self._ngood = 0
        self._nbad  = 0
    
    def close_enough(self, site):
        if not self.current_sites:
            return True
        if self.maxdist:
            d = int(site[1]) - int(self.current_sites[-1][2])
            # print(f"{self.current_sites[-1][2]} - {site[1]}: distance={d}")
            return d <= self.maxdist
        return True
    
    def add(self, site):
        # print(site)
        chrom = site[0]
        if chrom != self.current_chrom:
            if self.in_dmr:
                self.report_dmr()
            self.reset()
            self.current_chrom = site[0]
        dm = float(site[3])
        
        if self.in_dmr:
            if dm * self.direction > self.threshold:      # over threshold in wanted direction?
                if self.close_enough(site):               # if not too far away
                    # print(f"{site[1]} - extending DMR")
                    self.current_sites.append(site)       # store this site
                    self._ngood += 1                          # increase number of good sites
                    self._nbad = 0                            # and reset bad sites counter
                else:
                    # print(f"{site[1]} - too far, closing DMR")
                    self.report_dmr()
                    self.reset()
                    self.in_dmr = False
            else:
                self._nbad += 1
                if self._nbad > self.nskip:  # too many consecutive bad sites?
                    # print(f"{self._nbad} bad seen, closing DMR")
                    self.report_dmr()
                    self.reset()
                    self.in_dmr = False      # close DMR 
        else:
            if dm * self.direction > self.threshold:     # over threshold in wanted direction?
                if self.close_enough(site):
                    self.current_sites.append(site)      # store this site
                    self._ngood += 1             # and keep track of how many we saw
                    if self._ngood == self.nopen: # if we saw enough good sites...
                        # print(f"{self._ngood} seen, opening DMR")
                        self.in_dmr = True       # open DMR

    def report_dmr(self):
        #print(self.current_sites)
        if not self.current_sites:
            return
        start  = int(self.current_sites[0][1])
        end    = int(self.current_sites[-1][2])

        if self.slop > 1:
            start = start - self.slop
            end   = end   + self.slop
        elif self.slop > 0:
            p = int((end-start) * self.slop)
            start = start - p
            end   = end   + p

        size   = end - start
        nsites = len(self.current_sites)
        meth   = [s[3] for s in self.current_sites]
        mean   = stats.mean(meth)
        pstdev = stats.pstdev(meth)
        self.outstream.write(f"{self.current_chrom}\t{start}\t{end}\t{size}\t{nsites}\t{mean}\t{pstdev}\n")

    def close(self):
        if self.in_dmr:
            self.report_dmr()
        self.outstream.close()


class WindowCountDMRfinder(DMRfinder):
    """
    Count the number of hyper- and hypo-methylated sites (defined by threshold) in each window, and report
    windows in which the ratio of the two values is greater than ratio.
    """
    direction = 1
    winsize   = 10000
    ratio     = 2
    threshold = 0.2

    _winstart = 0
    _winend   = 0
    _nup      = 0
    _ndown    = 0

    def reset(self):
        self._nup = 0
        self._ndown = 0
        self._winend = self._winstart + self.winsize - 1
    
    def add(self, site):
        # print(site)
        chrom = site[0]
        if chrom != self.current_chrom:
            self.report_dmr()
            self.current_chrom = chrom
            self._winstart = 0
            self.reset()
            
        pos = int(site[1])
        dm  = site[3]
        if pos <= self._winend:              # still inside window?
            if dm > self.threshold:
                self._nup += 1
            elif dm < -self.threshold:
                self._ndown += 1
        else:
            # we're outside of window
            self.report_dmr()
            self._winstart = self.winsize * int(pos / self.winsize)
            self.reset()

    def report_dmr(self):
        #print(f"{self._winstart} {self._winend} {self._nup} {self._ndown}")
        
        def write(r):
            self.outstream.write(f"{self.current_chrom}\t{self._winstart}\t{self._winend}\t{self._nup}\t{self._ndown}\t{r}\n")
        
        if self._nup + self._ndown == 0:
            return
        
        if self.direction == 1:
            if self._ndown:
                r = self._nup / self._ndown
                if r >= self.ratio:
                    write(r)
            else:
                write("-")
        else:
            if self._nup:
                r = self._ndown / self._nup
                if r >= self.ratio:
                    write(r)
            else:
                write("-")

    def close(self):
        self.report_dmr()
        self.outstream.close()
            

DETECTOR_CLASSES = {
    'StreamingDMRfinder': StreamingDMRfinder,
    'WindowCountDMRfinder': WindowCountDMRfinder
}

class DMRmanager():
    detectors = []
    test_col  = None
    ctrl_col  = None
    diff_col  = None

    def init(self, configfile):
        sys.stderr.write(f'Reading config file {configfile}...\n')
        self.detectors = []
        d = None
        detname = ""
        with open(configfile, "r") as f:
            for line in f:
                if line.startswith('#'):
                    continue
                if line.startswith('--'):
                    det = None
                    continue
                words = [ w.strip() for w in line.split(':') ]
                if len(words) != 2:
                    continue
                key = words[0]
                if '_' in key:
                    tag = key[-1]
                    key = key[:-2]
                else:
                    tag = 's'
                if tag == 's':
                    value = words[1]
                elif tag == 'i':
                    value = int(words[1])
                elif tag == 'f':
                    value = float(words[1])
                elif tag == 'b':
                    value = (words[1] == 'True')
                    
                if key == 'type':
                    if value in DETECTOR_CLASSES:
                        detname = value
                        d = DETECTOR_CLASSES[value]() # initialize detector onbject
                        self.detectors.append(d)
                    else:
                        sys.stderr.write(f'Unknown detector type: {value}.\n')
                        sys.exit(2)
                else:
                    if hasattr(d, key):
                        setattr(d, key, value)
                    else:
                        sys.stderr.write(f'Detector {detname} does not have an attribute called {key}.\n')
                        sys.exit(3)
        sys.stderr.write(f'{len(self.detectors)} detectors defined.\n')
        
    def run(self, infile):
        for d in self.detectors:
            d.open()

        with open(infile, "rt") as f:
            c = csv.reader(f)
            next(c)
            for row in c:
                if self.diff_col:
                    dm = float(row[self.diff_col])
                else:
                    m_test = float(row[self.test_col])
                    m_ctrl = float(row[self.ctrl_col])
                    dm = (m_test - m_ctrl) / 100        

                site = [row[0], row[1], row[2], dm]
                #print(site)
                for d in self.detectors:
                    d.add(site)

            for d in self.detectors:
                d.close()

def main():
    parser = argparse.ArgumentParser(description="Detect DMRs from methylation data.",
                                     epilog="""See the comments in the provided configuration file example (config.txt) for a description of its format.\n\n
The input file should be comma-delimited (or tab-delimited if --tab is specified) and should have chromosome, start, and end as the first three columns.""",
                                    formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("configfile", help="Configuration file containing definition of DMR detectors.")
    parser.add_argument("infile", help="File containing methylation data.")
    parser.add_argument("--test_col", type=int, default=None, help="Column containing methylation values for test sample (1-based)")
    parser.add_argument("--ctrl_col", type=int, default=None, help="Column containing methylation values for control sample (1-based)")
    parser.add_argument("--diff_col", type=int, default=None, help="Column containing differential methylation values (1-based)")
    parser.add_argument("--tab", action="store_true", help="If specified, input file is tab-delimited instead of comma-delimited.")
    args = parser.parse_args()

    M = DMRmanager()
    if args.diff_col:
        M.diff_col = args.diff_col - 1
    elif args.test_col and args.ctrl_col:
        M.test_col = args.test_col - 1
        M.ctrl_col = args.ctrl_col - 1
    else:
        sys.stderr.write("Error: please specify either --diff_col or both --test_col and --ctrl_col.\n")
        sys.exit(1)

    M.init(args.configfile)
    M.run(args.infile)

if __name__ == "__main__":
    main()
    
