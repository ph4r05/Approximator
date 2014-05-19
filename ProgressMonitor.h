/* 
 * File:   ProgressMonitor.h
 * Author: ph4r05
 *
 * Created on May 19, 2014, 1:52 PM
 */

#ifndef PROGRESSMONITOR_H
#define	PROGRESSMONITOR_H

class ProgressMonitor {
private:    
    double step;
    double last;
    
public:
    ProgressMonitor(double step1) : step(step1), last(0.0) {};
    virtual ~ProgressMonitor();

    void setCur(double current);
    void reset() { last=0.0; };
};

#endif	/* PROGRESSMONITOR_H */

