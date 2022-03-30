###################################################
Timespan class & other utilities (`astropack.util`)
###################################################

**************************************
Getting Started with Class |Timespan|
**************************************

>>> from datetime import datetime, timezone, timedelta
>>> from astropy.time import Time, TimeDelta
>>> from astropack.util import Timespan
>>> dt1 = datetime(2016, 9, 10, 0, 0, 0, tzinfo=timezone.utc)
>>> dt2 = dt1 + timedelta(hours=1.5)
>>> ts1 = Timespan(Time(dt1), Time(dt2))  # input times must hve UTC timezone.
>>> ts1
Timespan(2016-09-10 00:00:00.000, 2016-09-10 01:30:00.000)
>>> ts1.start, ts1.end, ts1.midpoint
(<Time object: scale='utc' format='datetime' value=2016-09-10 00:00:00>,
 <Time object: scale='utc' format='datetime' value=2016-09-10 01:30:00>,
 <Time object: scale='utc' format='datetime' value=2016-09-10 00:45:00>)
>>> ts1.seconds, ts1.days
(5400.0, 0.0625)
>>> ts2 = ts1.delay_by(600)  # seconds
>>> ts2
Timespan(2016-09-10 00:10:00.000, 2016-09-10 01:40:00.000)
>>> ts2.expand_by(3600)  # seconds
Timespan(2016-09-09 23:10:00.000, 2016-09-10 02:40:00.000)
>>> ts2.intersection(ts1)
Timespan(2016-09-10 00:10:00.000, 2016-09-10 01:30:00.000)
>>> ts2.union(ts1)
Timespan(2016-09-10 00:00:00.000, 2016-09-10 01:40:00.000)
>>> ts2.subtract(ts1)
Timespan(2016-09-10 01:30:00.000, 2016-09-10 01:40:00.000)
>>> ts1.subtract(ts2)
Timespan(2016-09-10 00:00:00.000, 2016-09-10 00:10:00.000)
>>> ts2.split_at(Time('2016-09-10 00:30:00'))
[Timespan(2016-09-10 00:10:00.000, 2016-09-10 00:30:00.000),
 Timespan(2016-09-10 00:30:00.000, 2016-09-10 01:40:00.000)]
>>> ts1.periodic_events(ref_time=Time('2016-09-10 00:30:00'), \
                        period=TimeDelta(1200, format='sec'))
[<Time object: scale='utc' format='iso' value=2016-09-10 00:10:00.000>,
 <Time object: scale='utc' format='iso' value=2016-09-10 00:30:00.000>,
 <Time object: scale='utc' format='iso' value=2016-09-10 00:50:00.000>,
 <Time object: scale='utc' format='iso' value=2016-09-10 01:10:00.000>,
 <Time object: scale='utc' format='iso' value=2016-09-10 01:30:00.000>]
>>> str(ts1)
"Timespan '2016-09-10 00:00:00' to '2016-09-10 01:30:00' = 5400.000 seconds."

***************
Reference/API
***************

.. automodapi:: astropack.util
   :no-inheritance-diagram:
