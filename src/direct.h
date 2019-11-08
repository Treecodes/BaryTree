#ifndef H_DIRECT_H
#define H_DIRECT_H

void directSummation(double *source_x, double *source_y, double *source_z, double *source_charge, double *source_weight,
                double *target_x, double *target_y, double *target_z, double *target_charge,
                int number_of_sources, int number_of_targets, double *potential, double *total_potential, double kernel_parameter,
                char *kernel_name);

#endif /* H_DIRECT_H */
