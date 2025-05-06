import pandas as pd
import numpy as np
import math
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.special
import scipy.stats as stats

def generate_extra_plots(session: "pybmds.Session"):
    """
    Generate nested litter plots comparing observed vs estimated responder distributions
    for each dose level.

    Args:
        session (pybmds.Session): A fitted BMDS session with nested models.

    Returns:
        matplotlib.figure.Figure: The final plot figure.
    """

    rows, rows2 = [], []

    for j, model in enumerate(session.models):
        data = {"bmds_model_index": j,
               }
        data2 = {"bmds_model_index": j,
               }
        if model.has_results:
            res = model.results
            name = model.name()
            data.update(
                {"parameters": res.parameter_names,
                 "parameter_vals": res.parameters,
                 "model_name": name,
                }
            )
            data2.update(
                {"dose": res.litter.dose,
                 "est_prob": res.litter.estimated_probabilities,
                 "lsc": res.litter.lsc,
                 "litter_size": res.litter.litter_size,
                "scaled_residuals": res.litter.scaled_residuals,
                 "obs": res.litter.observed,
                 "model_name": name,
                }
            )
            rows.append(data)
            rows2.append(data2)

    
    litter = []
    for model in rows2:
        doses = model['dose']
        est_probs = model['est_prob']
        lsces = model['lsc']
        litter_sizes = model['litter_size']
        residuals = model['scaled_residuals']
        obs = model['obs']
        
        for i in range(len(doses)):
            litter.append({
                'bmds_model_index': model['bmds_model_index'],
                'model_name': model['model_name'],
                'dose': doses[i],
                'est_prob': est_probs[i],
                'lsc': lsces[i],
                'litter_size': litter_sizes[i],
                'obs': obs[i],
                'scaled_residuals': residuals[i]
            })
    litter = pd.DataFrame(litter)

    parameters = []
    for model in rows:
        param_dict = {'bmds_model_index': model['bmds_model_index']}
        for param, value in zip(model['parameters'], model['parameter_vals']):
            param_dict[param] = value
        parameters.append(param_dict)
    df_parameters = pd.DataFrame(parameters)

    final = pd.merge(litter, df_parameters, on='bmds_model_index', how='left')

    unique_doses = sorted(final['dose'].unique())
    dose_to_phi = {dose: f'phi{index + 1}' for index, dose in enumerate(unique_doses)}

    mname = ['Nested Logistic (lsc+ilc+)', 'Nested Logistic (lsc-ilc+)', 'NCTR (lsc+ilc+)', 'NCTR (lsc-ilc+)']

    def apply_phi(row, mname):
        if row['model_name'] in mname:
            phi_col = dose_to_phi.get(row['dose'])
            if phi_col and phi_col in row and row[phi_col] != 0:
                return row['est_prob'] * (1 / row[phi_col] - 1)
        return None

    def calc_Rd(row, mname):
        if row['model_name'] in mname:
            alpha = row['alpha']
            beta = row['beta']
            litter_size = row['litter_size']

            phi_column = dose_to_phi.get(row['dose'])

            if phi_column and phi_column in row and row[phi_column] == 0:
                return None 

            if any(x is None or x <= 0 or math.isnan(x) for x in [alpha, beta, litter_size]):
                return None 

            return (
                scipy.special.gamma(litter_size + 1) *
                scipy.special.gamma(alpha + beta) /
                (scipy.special.gamma(alpha) * scipy.special.gamma(beta) * scipy.special.gamma(litter_size + alpha + beta))
            )

        return None

    final['alpha'] = final.apply(lambda row: apply_phi(row, mname), axis=1)
    final['beta'] = final['alpha'] * (1 - final['est_prob']) / final['est_prob']
    final['Rd'] = final.apply(lambda row: calc_Rd(row, mname), axis=1)

    def generate_litter_columns(df, mname, dose_to_phi):
        """
        This function goes through the methodology described above for both ilc+ and ilc- models. If an ilc+ model has phi = 0, it will follow the ilc- approach.
        """
        max_litter_size = df['litter_size'].max()

        for i in range(int(max_litter_size) + 1):
            col_name = f"Litter_{i}" 
            
            def calc_first_value(row):
                phi_column = dose_to_phi.get(row['dose'], None)
                is_ilc_plus = row['model_name'] in mname
                phi_is_zero = phi_column and phi_column in row and row[phi_column] == 0

                if is_ilc_plus and not phi_is_zero:
                    return 1 
                else:
                    return 0 if i > int(row['litter_size']) else math.factorial(int(row['litter_size'])) / (math.factorial(i) * math.factorial(int(row['litter_size']) - i))
            
            df[col_name] = df.apply(calc_first_value, axis=1)

        for i in range(int(max_litter_size) + 1):
            col_name = f"Litter_Next_{i}"
            
            def calc_second_value(row):
                phi_column = dose_to_phi.get(row['dose'], None)
                is_ilc_plus = row['model_name'] in mname
                phi_is_zero = phi_column and phi_column in row and row[phi_column] == 0

                if is_ilc_plus and not phi_is_zero:
                    if i > row['litter_size']:
                        return 0
                    return scipy.special.gamma(i + row['alpha']) * scipy.special.gamma(row['litter_size'] - i + row['beta']) / (
                            scipy.special.gamma(i + 1) * scipy.special.gamma(row['litter_size'] - i + 1))
                else:
                    return (row['est_prob'] ** i) * ((1 - row['est_prob']) ** (row['litter_size'] - i))

            df[col_name] = df.apply(calc_second_value, axis=1)

        for i in range(int(max_litter_size) + 1):
            col_name = f"Litter_Final_{i}"

            def calc_final_value(row):
                phi_column = dose_to_phi.get(row['dose'], None)
                is_ilc_plus = row['model_name'] in mname
                phi_is_zero = phi_column and phi_column in row and row[phi_column] == 0

                if is_ilc_plus and not phi_is_zero:
                    return row[f"Litter_{i}"] * row[f"Litter_Next_{i}"] * row['Rd']
                else:
                    return row[f"Litter_{i}"] * row[f"Litter_Next_{i}"]

            df[col_name] = df.apply(calc_final_value, axis=1)
            
        return df

    df = generate_litter_columns(final, mname, dose_to_phi)

    def sum_litter_by_dose(df):
        """
        Estimated counts for each dose group
        """
        litter_final_cols = [col for col in df.columns if col.startswith("Litter_Final_")]
        grouped_df = df.groupby(['dose', 'model_name'])[litter_final_cols].sum().reset_index()

        return grouped_df

    def count_litter_obs(df):
        """
        Calculates the number of litters in each dose group that have exactly x responders
        """
        litter_final_cols = [col for col in df.columns if col.startswith("Litter_Final_")]

        result_df = pd.DataFrame()

        grouped_df = df.groupby(['dose', 'model_name'])

        result_df['dose'] = grouped_df['dose'].first().values
        result_df['model_name'] = grouped_df['model_name'].first().values

        for litter_col in litter_final_cols:
            litter_index = int(litter_col.split('_')[-1]) 
            
            counts = []

            for _, group in grouped_df:
                count = (group['obs'] == litter_index).sum()
                counts.append(count)

            result_df[litter_col] = counts

        result_df.reset_index(drop=True, inplace=True)

        return result_df

    def add_all_doses(counted_df):
        """
        Adds rows at the end for the estimated and observed values for all doses
        """
        summed_df = counted_df.groupby('model_name')[counted_df.columns.difference(['dose', 'model_name'])].sum()
        summed_df['dose'] = 'All'
        summed_df.reset_index(inplace=True)
        
        final_df = pd.concat([counted_df, summed_df], ignore_index=True)

        return final_df

    estimated = add_all_doses(sum_litter_by_dose(df))
    observed = add_all_doses(count_litter_obs(df))

    def transform(observed_df, estimated_df):
        obs_long = observed_df.melt(id_vars=['dose', 'model_name'], var_name='Litter_Final', value_name='obs_count')
        obs_long['Litter_Final'] = obs_long['Litter_Final'].str.extract(r'(\d+)').astype(int)

        est_long = estimated_df.melt(id_vars=['dose', 'model_name'], var_name='Litter_Final', value_name='est_count')
        est_long['Litter_Final'] = est_long['Litter_Final'].str.extract(r'(\d+)').astype(int)

        return pd.merge(obs_long, est_long, on=['dose', 'model_name', 'Litter_Final'], how='outer')

    plot_data = transform(observed, estimated)

    def calculate_plotting(incidence, n, alpha=0.05):
        p = incidence / float(n)
        z = stats.norm.ppf(1 - alpha / 2)
        z2 = z * z
        q = 1.0 - p
        tmp1 = 2 * n * p + z2
        ll = ((tmp1 - 1) - z * np.sqrt(z2 - (2 + 1 / n) + 4 * p * (n * q + 1))) / (2 * (n + z2))
        ul = ((tmp1 + 1) + z * np.sqrt(z2 + (2 + 1 / n) + 4 * p * (n * q - 1))) / (2 * (n + z2))
        return p, ll, ul

    plot_data['total_litters_per_dose'] = plot_data.groupby(['dose', 'model_name'])['obs_count'].transform('sum')
    plot_data[['obs_proportion', 'obs_ci_lower', 'obs_ci_upper']] = plot_data.apply(
        lambda row: calculate_plotting(row['obs_count'], row['total_litters_per_dose']), axis=1, result_type='expand'
    )
    plot_data['obs_ci_lower_count'] = plot_data['obs_ci_lower'] * plot_data['total_litters_per_dose']
    plot_data['obs_ci_upper_count'] = plot_data['obs_ci_upper'] * plot_data['total_litters_per_dose']

    sns.set(style="whitegrid")
    g = sns.FacetGrid(plot_data, col="dose", hue="model_name", height=5, aspect=2, col_wrap=1, sharey=False, sharex=False)
    g.map(sns.scatterplot, 'Litter_Final', 'obs_count', marker="o", label="Observed", s=50, color="black")

    for ax, dose in zip(g.axes.flat, plot_data['dose'].unique()):
        dose_data = plot_data[plot_data['dose'] == dose]
        means = dose_data['obs_count']
        lls = dose_data['obs_ci_lower_count']
        uls = dose_data['obs_ci_upper_count']
        ax.errorbar(dose_data['Litter_Final'], means,
                    yerr=[(means - lls).clip(0), (uls - means).clip(0)],
                    fmt='o', color='black', capsize=4)

    g.map(sns.lineplot, 'Litter_Final', 'est_count', marker="o", label="Estimated")
    g.add_legend(title="Model", label_order=["Observed"] + list(g._legend_data.keys()))
    g.set_axis_labels("Number of Responders", "Number of Litters")
    g.set_titles("Dose = {col_name}")
    for ax in g.axes.flat:
        ax.set_xlabel('Number of Responders')
        ax.tick_params(axis='x', rotation=45)
    plt.subplots_adjust(hspace=0.4)

    return g.fig

# import pybmds
# import pandas as pd
# import math

# df = pd.read_csv("Nested_datasets.csv")
# df = df[df['ID'] == 3]

# for id, rows in df.groupby('ID'):
#     dataset = pybmds.NestedDichotomousDataset(
#         name = "Chlorothalonil Implantations",
#         doses = rows.Dose.tolist(),
#         litter_ns = rows.N.tolist(),
#         incidences = rows.Incidence.tolist(),
#         litter_covariates = rows.lsc.tolist()
#     )

# session = pybmds.Session(dataset=dataset)
# session.add_model(pybmds.Models.NCTR, settings = {"litter_specific_covariate": 0, "intralitter_correlation": 0})
# #session.add_default_models(settings = {"litter_specific_covariate": 0, "intralitter_correlation": 0})
# #session.add_default_models(settings = {"litter_specific_covariate": 1, "intralitter_correlation": 1})
# #session.add_default_models(settings = {"litter_specific_covariate": 0, "intralitter_correlation": 1})
# #session.add_default_models(settings = {"litter_specific_covariate": 1, "intralitter_correlation": 0})
# session.execute()

# session.plot()

# fig = generate_extra_plots(session)

# fig.savefig('plots.png')