import matplotlib.pyplot as plt

def hstat(data, ax=None, **kwargs):
    if ax is None:
        fig, ax = plt.subplots(1, 1)
    else:
        fig = plt.gcf()

    avg = data.filter(data.is_not_nan()).mean()
    std = data.filter(data.is_not_nan()).std()

    if 'label' in kwargs:
        label = kwargs.pop('label')
    else:
        label = ""

    label = rf"{label} ($\mu={avg:.2f}$ ; $\sigma={std:.2f}$)"

    ax.hist(data, label=label, **kwargs)

    return fig