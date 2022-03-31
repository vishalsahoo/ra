// Simple static time logger class. Only takes in doubles. Not thread safe
// Contributed by Hans
class TimeLog
{
public:
    // Ignores the first num_ignored time variables entered. Must be called before using log
    static void ignore(int num_ignored)
    {
        TimeLog::num_ignored = num_ignored;
    }
    // Adds a time to the logger
    static void log(const double& time)
    {
        ++count;
        if (count > num_ignored)
            accum += time;
    }
    // Returns the average of the results
    static double calc_avg()
    {
        if (count > num_ignored)
            return (accum / (double)(count - num_ignored));
        else
            return -1;
    }
    // Resets the logger to its initial state
    static void reset()
    {
        num_ignored = 0;
        accum = 0;
        count = 0;
    }
protected:
    static double accum;
    static int num_ignored;
    static int count;
};
double TimeLog::accum = 0;
int TimeLog::num_ignored = 0;
int TimeLog::count = 0;