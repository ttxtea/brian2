{# IS_OPENMP_COMPATIBLE #}
{% macro cpp_file() %}

#include "network.h"
#include<stdlib.h>
#include<iostream>
#include <ctime>
#include<utility>
{{ openmp_pragma('include') }}

#define Clock_epsilon 1e-14

Network::Network()
{
	t = 0.0;
}

void Network::clear()
{
	objects.clear();
}

void Network::add(Clock* clock, codeobj_func func)
{
#if defined(_MSC_VER) && (_MSC_VER>=1700)
	objects.push_back(std::make_pair<Clock*, codeobj_func>(std::move(clock), std::move(func)));
#else
	objects.push_back(std::make_pair<Clock*, codeobj_func>(clock, func));
#endif
}

void Network::run(const double duration, void (*report_func)(const double, const double, const double), const double report_period)
{
    {% if openmp_pragma('with_openmp') %}
    double start;
    {% else %}
	std::clock_t start, current;
	{% endif %}
    const double t_start = t;
	const double t_end = t + duration;
	double next_report_time = report_period;
	// compute the set of clocks
	compute_clocks();
	// set interval for all clocks

	for(std::set<Clock*>::iterator i=clocks.begin(); i!=clocks.end(); i++)
		(*i)->set_interval(t, t_end);

    {% if openmp_pragma('with_openmp') %}
    start = omp_get_wtime();
    {% else %}
	start = std::clock();
    {% endif %}
	if (report_func)
	{
	    report_func(0.0, 0.0, duration);
	}

	Clock* clock = next_clocks();

	{{ openmp_pragma('parallel') }}
	{
		while(clock->running())
		{
			for(int i=0; i<objects.size(); i++)
			{
				{{ openmp_pragma('single-nowait') }}
				{
					if (report_func)
		            {
                        {% if openmp_pragma('with_openmp') %}
                        const double elapsed = omp_get_wtime() - start;
                        {% else %}
		                current = std::clock();
                        const double elapsed = (double)((current - start) / CLOCKS_PER_SEC);
                        {% endif %}
		                if (elapsed > next_report_time)
		                {
		                    report_func(elapsed, (clock->t_()-t_start)/duration, duration);
		                    next_report_time += report_period;
		                }
		            }
		        }
				Clock *obj_clock = objects[i].first;
				// Only execute the object if it uses the right clock for this step
				if (curclocks.find(obj_clock) != curclocks.end())
				{
	                codeobj_func func = objects[i].second;
	                {{ openmp_pragma('barrier') }}
	                func();
				}
			}
			{{ openmp_pragma('single') }}
			{
				for(std::set<Clock*>::iterator i=curclocks.begin(); i!=curclocks.end(); i++)
					(*i)->tick();
			}
			clock = next_clocks();
		}
		{{ openmp_pragma('single-nowait') }}
		{
			if (report_func)
			{
			    {% if openmp_pragma('with_openmp') %}
                report_func(omp_get_wtime() - start, 1.0, duration);
                {% else %}
			    current = std::clock();
			    report_func((double)(current - start)/({{ openmp_pragma('get_num_threads') }} * CLOCKS_PER_SEC), 1.0, duration);
			    {% endif %}
			}
		}
	t = t_end;
	}
}

void Network::compute_clocks()
{
	clocks.clear();
	for(int i=0; i<objects.size(); i++)
	{
		Clock *clock = objects[i].first;
		clocks.insert(clock);
	}
}

Clock* Network::next_clocks()
{
	// find minclock, clock with smallest t value
	Clock *minclock = *clocks.begin();
	for(std::set<Clock*>::iterator i=clocks.begin(); i!=clocks.end(); i++)
	{
		Clock *clock = *i;
		if(clock->t_()<minclock->t_())
			minclock = clock;
	}
	// find set of equal clocks
	{{ openmp_pragma('single') }}
	{
		// find set of equal clocks
		curclocks.clear();

		double t = minclock->t_();
		for(std::set<Clock*>::iterator i=clocks.begin(); i!=clocks.end(); i++)
		{
			Clock *clock = *i;
			double s = clock->t_();
			if(s==t || fabs(s-t)<=Clock_epsilon)
				curclocks.insert(clock);
		}
	}
	return minclock;
}

{% endmacro %}

{% macro h_file() %}

#ifndef _BRIAN_NETWORK_H
#define _BRIAN_NETWORK_H

#include<vector>
#include<utility>
#include<set>
#include "brianlib/clocks.h"

typedef void (*codeobj_func)();

class Network
{
	std::set<Clock*> clocks, curclocks;
	void compute_clocks();
	Clock* next_clocks();
public:
	std::vector< std::pair< Clock*, codeobj_func > > objects;
	double t;

	Network();
	void clear();
	void add(Clock *clock, codeobj_func func);
	void run(const double duration, void (*report_func)(const double, const double, const double), const double report_period);
};

#endif

{% endmacro %}
