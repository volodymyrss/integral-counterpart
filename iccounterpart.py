import counterpart

class Event(counterpart.da.DataAnalysis):
    run_for_hashe=True

    event_kind=counterpart.IceCubeEvent

    def main(self):
        return self.event_kind