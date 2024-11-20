#!/usr/bin/env python3


class SignalCollection:
    """
    Collects processed dataframes (signals), types of data (sigtypes:
    'continuous' or 'binary'), along with relevant lists (lists: e.g.
    intervals between budding events).  Types of data informs the processes
    that are applied to them.
    """

    def __init__(self) -> None:
        self.signals = {}
        # 'continous' or 'binary'
        self.sigtypes = {}
        self.lists = {}
