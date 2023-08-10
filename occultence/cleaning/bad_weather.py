

def mask_bad_weather(self):

    print(str(100 * np.divide(float(np.count_nonzero(self.mask_bad_weather)),
                              len(self.mask_bad_weather))) + "% of data is flagged as bad weather")

