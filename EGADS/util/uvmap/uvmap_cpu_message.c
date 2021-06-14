#include "UVMAP_LIB.h"

/*
 * UVMAP : TRIA-FACE SURFACE MESH UV MAPPING GENERATOR
 *         DERIVED FROM AFLR4, UG, UG2, and UG3 LIBRARIES
 * $Id: uvmap_cpu_message.c,v 1.3 2020/07/04 20:46:56 marcum Exp $
 * Copyright 1994-2020, David L. Marcum
 */

#if defined (__APPLE__) && defined (__MACH__)
#include <mach/clock.h>
#include <mach/mach.h>
#endif

static double uvmap_clock_gettime (void);

typedef struct _cpu_timer_struct cpu_timer_struct;

struct _cpu_timer_struct 
{
  char CPU_Timer_Label[512];
  double CPU_Time_Start;
  double CPU_Time_Total;
  double CPU_Time_Total_Overall;
};

void uvmap_cpu_message (char *Label) {

  // Write CPU time since last call to standard output.

  char Text[512], Text_Label[512];

  static double initial_CPU_Time = -1;
  static double old_CPU_Time = -1;

  double CPU_Time, new_CPU_Time;

  if (strcmp (Label, "reset") == 0) {
    old_CPU_Time = -1;
    return;
  }

  new_CPU_Time = uvmap_clock_gettime ();

  if (old_CPU_Time < 0)
    old_CPU_Time = new_CPU_Time;

  CPU_Time = new_CPU_Time - old_CPU_Time;
  CPU_Time = MAX (CPU_Time, 0);

  old_CPU_Time = new_CPU_Time;

  if (strcmp (Label, "start_overall") == 0)
  {
    initial_CPU_Time = new_CPU_Time;

    old_CPU_Time = -1.0;

    return;
  }

  else if (strcmp (Label, "") == 0)
    return;

  else if (strcmp (Label, "end_overall") == 0)
  {
    CPU_Time = new_CPU_Time - initial_CPU_Time;

    strcpy (Text_Label, "OVERALL  :");

    initial_CPU_Time = -1.0;

    old_CPU_Time = -1.0;
  }
  else
    strcpy (Text_Label, Label);

  if (CPU_Time < 60)
    snprintf (Text, 512, "%s CPU Time          =%10.3f   seconds",
              Text_Label, CPU_Time);
  else if (CPU_Time < 3600)
    snprintf (Text, 512, "%s CPU Time          =%10.3f   minutes",
              Text_Label, CPU_Time/60.0);
  else
    snprintf (Text, 512, "%s CPU Time          =%10.3f     hours",
              Text_Label, CPU_Time/3600.0);

  uvmap_message (Text);

  return;
}

void uvmap_cpu_timer (char *Task, char * Label) {

  // Time multiple calls of a labeled operation and output the total CPU time.

  char Text[512];
  char Text_Label[21];

  static cpu_timer_struct *cpu_timer_struct_ptr = NULL;

  cpu_timer_struct *cpu_timer_struct_ptr_i;

  static INT_ n_cpu_timers = 0;
  static INT_ cpu_timer_err = 0;

  INT_ i = 0;
  INT_ found = 0;

  double CPU_Time_Total;

  if (cpu_timer_err)
    return;

  if (strcmp (Task, "create") == 0)
  {
    while (i < n_cpu_timers && found == 0)
    {
      cpu_timer_struct_ptr_i = &(cpu_timer_struct_ptr[i]);

      if (strcmp (cpu_timer_struct_ptr_i->CPU_Timer_Label, Label) == 0)
      {
        found = 1;

        cpu_timer_struct_ptr_i->CPU_Time_Start = uvmap_clock_gettime();
      }

      ++i;
    }

    if (found == 0)
    {
      cpu_timer_struct_ptr = (cpu_timer_struct *) uvmap_realloc (&cpu_timer_err,
                                                              cpu_timer_struct_ptr,
                                                              (n_cpu_timers+1) * (sizeof (cpu_timer_struct)));

      if (cpu_timer_err)
        return;

      cpu_timer_struct_ptr_i = &(cpu_timer_struct_ptr[n_cpu_timers]);

      strcpy (cpu_timer_struct_ptr_i->CPU_Timer_Label, Label);

      cpu_timer_struct_ptr_i->CPU_Time_Total = 0;
      cpu_timer_struct_ptr_i->CPU_Time_Total_Overall = 0;

      ++n_cpu_timers;
    }
  }

  if (strcmp (Task, "cleanup") == 0 || cpu_timer_err)
  {
    uvmap_free (cpu_timer_struct_ptr);

    n_cpu_timers = 0;

    cpu_timer_struct_ptr = NULL;
  }

  else if (strcmp (Task, "start") == 0)
  {
    while (i < n_cpu_timers && found == 0)
    {
      cpu_timer_struct_ptr_i = &(cpu_timer_struct_ptr[i]);

      if (strcmp (cpu_timer_struct_ptr_i->CPU_Timer_Label, Label) == 0)
      {
        found = 1;

        cpu_timer_struct_ptr_i->CPU_Time_Start = uvmap_clock_gettime();
      }

      ++i;
    }
  }

  else if (strcmp (Task, "stop") == 0)
  {
    while (i < n_cpu_timers && found == 0)
    {
      cpu_timer_struct_ptr_i = &(cpu_timer_struct_ptr[i]);

      if (strcmp (cpu_timer_struct_ptr_i->CPU_Timer_Label, Label) == 0)
      {
        found = 1;

        cpu_timer_struct_ptr_i->CPU_Time_Total = cpu_timer_struct_ptr_i->CPU_Time_Total + uvmap_clock_gettime() - cpu_timer_struct_ptr_i->CPU_Time_Start;
        cpu_timer_struct_ptr_i->CPU_Time_Total_Overall = cpu_timer_struct_ptr_i->CPU_Time_Total_Overall + uvmap_clock_gettime() - cpu_timer_struct_ptr_i->CPU_Time_Start;
      }

      ++i;
    }
  }

  else if (strcmp (Task, "output") == 0)
  {
    if (strcmp (Label, "ALL") == 0)
      found = 2;
    if (strcmp (Label, "OVERALL") == 0)
      found = 3;

    while (i < n_cpu_timers && found != 1)
    {
      cpu_timer_struct_ptr_i = &(cpu_timer_struct_ptr[i]);

      if (found == 0 && strcmp (Label, cpu_timer_struct_ptr_i->CPU_Timer_Label) == 0)
        found = 1;

      if (found)
      {
        strcpy (Text_Label, "");
        strncat (Text_Label, cpu_timer_struct_ptr_i->CPU_Timer_Label, 20); 

        if (found <= 2)
          CPU_Time_Total = cpu_timer_struct_ptr_i->CPU_Time_Total;
        else
          CPU_Time_Total = cpu_timer_struct_ptr_i->CPU_Time_Total_Overall;

        if (CPU_Time_Total < 60)
          snprintf (Text, 512, "%-20s CPU Time=%10.3f   seconds",
                    Text_Label, CPU_Time_Total);
        else if (CPU_Time_Total < 3600)
          snprintf (Text, 512, "%-20s CPU Time=%10.3f   minutes",
                    Text_Label, CPU_Time_Total/60.0);
        else
          snprintf (Text, 512, "%-20s CPU Time=%10.3f     hours",
                    Text_Label, CPU_Time_Total/3600.0);

        uvmap_message (Text);

        if (found)
          cpu_timer_struct_ptr_i->CPU_Time_Total = 0;

        if (found == 3)
          cpu_timer_struct_ptr_i->CPU_Time_Total_Overall = 0;
      }

      ++i;
    }
  }

  return;

}

double uvmap_clock_gettime (void) {

  // Get CPU clock time.

  double CPU_Time = 0;

#if defined (__APPLE__) && defined (__MACH__)

  clock_serv_t cclock;
  mach_timespec_t ts;

  host_get_clock_service (mach_host_self(), SYSTEM_CLOCK, &cclock);
  clock_get_time (cclock, &ts);
  mach_port_deallocate (mach_task_self(), cclock);

  CPU_Time = (double) (ts.tv_sec + ts.tv_nsec*1e-9);

#elif defined (__linux__)

  struct timespec ts;

  clock_gettime (CLOCK_MONOTONIC, &ts);

  CPU_Time = (double) (ts.tv_sec + ts.tv_nsec*1e-9);

#else

  CPU_Time = (double) clock () / ((double) CLOCKS_PER_SEC);

#endif

  return CPU_Time;
}
