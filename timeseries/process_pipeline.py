#!/usr/bin/env python3


def apply_postprocesses(signalcollection, process_dict):
    """Apply postprocesses

    Parameters
    ----------
    signalcollection : signalcollection.SignalCollection object
        Collect of signals to apply processes to
    process_dict : dict
        Dictionary of postprocesses to apply. Keys: string to append to name of
        source signal. Values: dictionary that specifies process, parameters,
        type of signal to operate on, type of signal of output.

    Examples
    --------
    FIXME: Add docs.

    """
    # Goes through each appending string, which represents a 'key' for each post-
    # process
    for process_appending_string, options_dict in process_dict.items():
        # Defines a list of signal (names) that I want my post-process to run on,
        # based on the ending of the signal name.
        applied_signames = [
            signame
            for signame in list(signalcollection.signals.keys())
            if signame.endswith(options_dict["signame_endswith"])
        ]
        for signame in applied_signames:
            # New signal name: append string that indicates which post-process
            # created it
            newname = signame + process_appending_string
            # Run if the signal type (continuous or binary) is correct
            if signalcollection.sigtypes[signame] == options_dict["input_sigtype"]:
                signalcollection.signals[newname] = options_dict["runner"].run(
                    signalcollection.signals[signame]
                )
                # Specify signal type of new signal
                signalcollection.sigtypes[newname] = options_dict["output_sigtype"]
            else:
                print(
                    f"{signame} is not {options_dict['input_sigtype']}, skipping {process_appending_string}"
                )
