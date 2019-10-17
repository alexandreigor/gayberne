#ifndef RAN1_H_
#define RAN1_H_

float ran0(long *idum);
float ran1(long *idum);

int getNTab();
void getRan1Status(long *_iy, long _iv[]);
void reinitRan1(const long _iy, const long _iv[]);

#endif /* RAN1_H_ */
