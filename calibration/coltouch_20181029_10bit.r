  0
 47
 64
 76
 87
 96
104
111
118
124
130
136
141
146
151
156
160
165
169
173
177
181
185
189
192
196
199
203
206
209
212
215
219
222
225
227
230
233
236
239
241
244
247
249
252
254
257
259
262
264
267
269
271
274
276
278
280
283
285
287
289
291
293
296
298
300
302
304
306
308
310
312
314
316
318
319
321
323
325
327
329
331
332
334
336
338
339
341
343
345
346
348
350
352
353
355
357
358
360
361
363
365
366
368
369
371
373
374
376
377
379
380
382
383
385
386
388
389
391
392
394
395
397
398
400
401
402
404
405
407
408
409
411
412
414
415
416
418
419
420
422
423
424
426
427
428
430
431
432
434
435
436
438
439
440
441
443
444
445
446
448
449
450
451
453
454
455
456
457
459
460
461
462
464
465
466
467
468
469
471
472
473
474
475
476
478
479
480
481
482
483
484
486
487
488
489
490
491
492
493
494
496
497
498
499
500
501
502
503
504
505
506
507
509
510
511
512
513
514
515
516
517
518
519
520
521
522
523
524
525
526
527
528
529
530
531
532
533
534
535
536
537
538
539
540
541
542
543
544
545
546
547
548
549
550
551
552
553
554
555
556
557
558
559
560
560
561
562
563
564
565
566
567
568
569
570
571
572
573
573
574
575
576
577
578
579
580
581
582
583
583
584
585
586
587
588
589
590
591
591
592
593
594
595
596
597
598
598
599
600
601
602
603
604
604
605
606
607
608
609
610
610
611
612
613
614
615
615
616
617
618
619
620
621
621
622
623
624
625
625
626
627
628
629
630
630
631
632
633
634
634
635
636
637
638
638
639
640
641
642
642
643
644
645
646
646
647
648
649
650
650
651
652
653
653
654
655
656
657
657
658
659
660
660
661
662
663
663
664
665
666
667
667
668
669
670
670
671
672
673
673
674
675
676
676
677
678
679
679
680
681
681
682
683
684
684
685
686
687
687
688
689
690
690
691
692
692
693
694
695
695
696
697
698
698
699
700
700
701
702
703
703
704
705
705
706
707
707
708
709
710
710
711
712
712
713
714
714
715
716
717
717
718
719
719
720
721
721
722
723
723
724
725
726
726
727
728
728
729
730
730
731
732
732
733
734
734
735
736
736
737
738
738
739
740
740
741
742
742
743
744
744
745
746
746
747
748
748
749
750
750
751
752
752
753
754
754
755
755
756
757
757
758
759
759
760
761
761
762
763
763
764
765
765
766
766
767
768
768
769
770
770
771
772
772
773
773
774
775
775
776
777
777
778
778
779
780
780
781
782
782
783
783
784
785
785
786
787
787
788
788
789
790
790
791
791
792
793
793
794
794
795
796
796
797
798
798
799
799
800
801
801
802
802
803
804
804
805
805
806
807
807
808
808
809
810
810
811
811
812
813
813
814
814
815
816
816
817
817
818
818
819
820
820
821
821
822
823
823
824
824
825
825
826
827
827
828
828
829
830
830
831
831
832
832
833
834
834
835
835
836
836
837
838
838
839
839
840
840
841
842
842
843
843
844
844
845
846
846
847
847
848
848
849
849
850
851
851
852
852
853
853
854
855
855
856
856
857
857
858
858
859
860
860
861
861
862
862
863
863
864
864
865
866
866
867
867
868
868
869
869
870
870
871
872
872
873
873
874
874
875
875
876
876
877
878
878
879
879
880
880
881
881
882
882
883
883
884
885
885
886
886
887
887
888
888
889
889
890
890
891
891
892
892
893
894
894
895
895
896
896
897
897
898
898
899
899
900
900
901
901
902
902
903
903
904
905
905
906
906
907
907
908
908
909
909
910
910
911
911
912
912
913
913
914
914
915
915
916
916
917
917
918
918
919
919
920
920
921
921
922
922
923
923
924
924
925
925
926
927
927
928
928
929
929
930
930
931
931
932
932
933
933
934
934
935
935
936
936
937
937
938
938
939
939
939
940
940
941
941
942
942
943
943
944
944
945
945
946
946
947
947
948
948
949
949
950
950
951
951
952
952
953
953
954
954
955
955
956
956
957
957
958
958
959
959
960
960
960
961
961
962
962
963
963
964
964
965
965
966
966
967
967
968
968
969
969
970
970
971
971
972
972
972
973
973
974
974
975
975
976
976
977
977
978
978
979
979
980
980
980
981
981
982
982
983
983
984
984
985
985
986
986
987
987
987
988
988
989
989
990
990
991
991
992
992
993
993
994
994
994
995
995
996
996
997
997
998
998
999
999
1000
1000
1000
1001
1001
1002
1002
1003
1003
1004
1004
1005
1005
1005
1006
1006
1007
1007
1008
1008
1009
1009
1010
1010
1010
1011
1011
1012
1012
1013
1013
1014
1014
1014
1015
1015
1016
1016
1017
1017
1018
1018
1019
1019
1019
1020
1020
1021
1021
1022
1022
1023
1023