
def debug_output():
    import dataanalysis.printhook
    dataanalysis.printhook.global_all_output=True
    dataanalysis.printhook.global_permissive_output=True
    dataanalysis.printhook.global_fancy_output=True
    dataanalysis.printhook.LogStreams=[]



def test_ligomap():
    import counterpart as ct

    le=ct.LIGOEvent(
            use_trigger_time="2017-06-08T02:01:16.492188",
            use_gname="G288732",
            use_loc_map_path="/home/savchenk/work/ligo/lvt170608/bayestar.fits"
        )

    le.get()

    a=le.loc_map

    assert abs(le.loc_region.max()-100)<1e-5


def test_visibility():
    import dataanalysis
    debug_output()

    import counterpart as ct

    le = ct.LIGOEvent(
        use_trigger_time="2017-06-08T02:01:16.492188",
        use_gname="G288732",
        use_loc_map_path="/home/savchenk/work/ligo/lvt170608/bayestar.fits"
    )
    le.promote()

    iv=ct.INTEGRALVisibility(use_minsolarangle=40).get()
    assert iv.minsolarangle == 40
    assert abs(iv.total_visible - 0.244737004674) < 1e-5

    iv = ct.INTEGRALVisibility(use_minsolarangle=50).get(update=True)

    assert  iv.minsolarangle == 50
    assert abs(iv.total_visible - 0.566222565777) < 1e-5

    assert iv.peak_target_visible() == (110.56640625, 9.3671571526015782)

def test_count_limits():
    import dataanalysis
    debug_output()

    import counterpart as ct

    le = ct.LIGOEvent(
        use_trigger_time="2017-06-08T02:01:16.492188",
        use_gname="G288732",
        use_loc_map_path="/home/savchenk/work/ligo/lvt170608/bayestar.fits"
    ).get()

    cl=ct.CountLimits().get()
    print cl.count_limits

def test_responses():
    import dataanalysis
    debug_output()

    import counterpart as ct

    le = ct.LIGOEvent(
        use_trigger_time="2017-06-08T02:01:16.492188",
        use_gname="G288732",
        use_loc_map_path="/home/savchenk/work/ligo/lvt170608/bayestar.fits"
    ).get()

    cl=ct.Responses().get()
    print cl.responses
