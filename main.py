# -*- coding: utf-8 -*-
import math
import time
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Literal

import cv2
import numpy as np

from SFR import SFR

ROI_SIZE = 70
AOI_LENGTH = 30
AOI_WIDTH = 20

orig_img = cv2.imread((Path('edge_samples') / '0007.png').as_posix())
start_time = time.time()
is_gray_image = True if len(orig_img.shape) == 2 else False

img = np.copy(orig_img)
# if is_gray_image:
#     img = np.copy(orig_img)
# else:
#     img = cv2.cvtColor(orig_img, cv2.COLOR_BGR2GRAY)

# 图片尺寸
width = img.shape[1]
height = img.shape[0]
diagonal_length = math.sqrt(width**2 + height**2)

# 计算中心坐标
center_x = width // 2
center_y = height // 2

cos_height = (height**2 + diagonal_length**2 - width**2) / (
    2 * height * diagonal_length
)
cos_width = (width**2 + diagonal_length**2 - height**2) / (2 * width * diagonal_length)


def distance(x1, y1, x2, y2):
    return np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)


@dataclass
class AOI:
    '''储存ROI中的T或S方向数据'''

    lt_pos: tuple[int, int]
    rt_pos: tuple[int, int]
    img: cv2.typing.MatLike


@dataclass
class ROI:
    '''储存各视场各方向数据'''

    img: cv2.typing.MatLike
    lt_pos: tuple[int, int]
    rt_pos: tuple[int, int]
    direction: (
        Literal['C'] | Literal['UL'] | Literal['UR'] | Literal['BL'] | Literal['BR']
    )
    field: float
    t_aoi: AOI | None = None
    s_aoi: AOI | None = None


# 计算中心视场需要计算MTF的区域
rois: list[ROI] = [
    ROI(
        img=orig_img[
            center_y - ROI_SIZE : center_y + ROI_SIZE,
            center_x - ROI_SIZE : center_x + ROI_SIZE,
        ],
        lt_pos=(center_x - ROI_SIZE, center_y - ROI_SIZE),
        rt_pos=(center_x + ROI_SIZE, center_y + ROI_SIZE),
        direction='C',
        field=0,
    )
]

# 获得每个视场每个方向需要计算MTF的区域
for field in (0.3, 0.5, 0.75, 0.9):
    # 计算偏移向量
    offset = field * (diagonal_length / 2)
    offset_x = int(offset * cos_width)
    offset_y = int(offset * cos_height)

    # 图像中间左上角坐标
    pos_ul = 'UL', (center_x - offset_x, center_y - offset_y)
    pos_ur = 'UR', (center_x + offset_x, center_y - offset_y)
    pos_bl = 'BL', (center_x - offset_x, center_y + offset_y)
    pos_br = 'BR', (center_x + offset_x, center_y + offset_y)

    for _direction, (_lt, _rb) in (pos_ul, pos_ur, pos_bl, pos_br):
        rois.append(
            ROI(
                img=orig_img[
                    _rb - ROI_SIZE : _rb + ROI_SIZE, _lt - ROI_SIZE : _lt + ROI_SIZE
                ],
                lt_pos=(_lt - ROI_SIZE, _rb - ROI_SIZE),
                rt_pos=(_lt + ROI_SIZE, _rb + ROI_SIZE),
                direction=_direction,
                field=field,
            )
        )

rect_size = 0
delta = 0

# 遍历ROI找角点，并以此找到AOI（单个TS方向）
for i in rois:

    '''
    shi-Tomas算法是对Harris角点检测算法的改进,一般会比Harris算法得到更好的角点.
    Harris算法的角点响应函数是将矩阵M的行列式值与M的迹相减,利用差值判断是否为角点,
    后来shi和Tomas提出改进的方法是,若矩阵M的两个特征值中较小的一个大于阈值,则认为
    他是角点,即:
        R=min(a1,a2)

    API:
        corners=cv.goodFeaturesToTrack(image, maxcorners, qualityLevel, minDistance)
        参数:
            image: 输入的灰度图像
            maxCorners: 获取角点数的数目
            qualityLevel: 该参数指出最低可接受的角点质量水平,在0~1之间
            minDistance: 角点之间的最小欧氏距离,避免得到相邻特征点
        返回:
            corners:搜索到的角点,在这里所有低于质量水平的角点被排除,然后把合格的角点按照质量排序,
            然后将质量较好的角点附近(小于最小欧氏距离)的角点删除,最后找到maxCorners个角点返回

    '''
    corners = cv2.goodFeaturesToTrack(
        cv2.cvtColor(i.img, cv2.COLOR_BGR2GRAY), 4, 0.5, 5
    )
    # TODO: 占用时间长

    if i.direction == 'C':
        x1, y1 = corners[0].ravel()
        x2, y2 = corners[1].ravel()
        x3, y3 = corners[2].ravel()
        x4, y4 = corners[3].ravel()

        # test=========================
        cv2.circle(
            img, (i.lt_pos[0] + int(x1), i.lt_pos[1] + int(y1)), 2, (0, 0, 255), 2
        )
        cv2.circle(
            img, (i.lt_pos[0] + int(x2), i.lt_pos[1] + int(y2)), 2, (0, 255, 0), 2
        )
        cv2.circle(
            img, (i.lt_pos[0] + int(x3), i.lt_pos[1] + int(y3)), 2, (255, 0, 0), 2
        )
        cv2.circle(
            img, (i.lt_pos[0] + int(x4), i.lt_pos[1] + int(y4)), 2, (255, 255, 0), 2
        )
        # test=========================

        distances: list[float] = [
            distance(x1, y1, x2, y2),
            distance(x1, y1, x3, y3),
            distance(x1, y1, x4, y4),
            distance(x2, y2, x3, y3),
            distance(x2, y2, x4, y4),
            distance(x3, y3, x4, y4),
        ]
        # 棋盘大小
        rect_size = int(np.min(distances))

        AOI_LENGTH = int(rect_size / 3 * 2)
        AOI_WIDTH = int(rect_size / 5 * 2)

        delta = int(
            np.min(np.abs([y2 - y1, y3 - y1, y4 - y1, x2 - x1, x3 - x1, x4 - x1]))
        )

        print(delta)

    x, y = corners[0].ravel()
    x, y = int(x), int(y)

    # test=========================
    x1, y1 = corners[0].ravel()
    x2, y2 = corners[1].ravel()
    x3, y3 = corners[2].ravel()
    x4, y4 = corners[3].ravel()
    cv2.circle(img, (i.lt_pos[0] + int(x1), i.lt_pos[1] + int(y1)), 2, (0, 0, 255), 2)
    cv2.circle(img, (i.lt_pos[0] + int(x2), i.lt_pos[1] + int(y2)), 2, (0, 255, 0), 2)
    cv2.circle(img, (i.lt_pos[0] + int(x3), i.lt_pos[1] + int(y3)), 2, (255, 0, 0), 2)
    cv2.circle(img, (i.lt_pos[0] + int(x4), i.lt_pos[1] + int(y4)), 2, (255, 255, 0), 2)
    # test=========================

    # 防止ROI选出ROI识别区域
    # 竖直方向超出
    if i.lt_pos[1] + y - (rect_size / 2) - AOI_LENGTH / 2 > i.lt_pos[1]:
        _lt_pos = (
            int(i.lt_pos[0] + x - (delta / 2) - AOI_WIDTH / 2),
            int(i.lt_pos[1] + y - (rect_size / 2) - AOI_LENGTH / 2),
        )
        _rb_pos = (
            int(i.lt_pos[0] + x - (delta / 2) + AOI_WIDTH / 2),
            int(i.lt_pos[1] + y - (rect_size / 2) + AOI_LENGTH / 2),
        )
    else:
        _lt_pos = (
            int(i.lt_pos[0] + x + (delta / 2) - AOI_WIDTH / 2),
            int(i.lt_pos[1] + y + (rect_size / 2) - AOI_LENGTH / 2),
        )
        _rb_pos = (
            int(i.lt_pos[0] + x + (delta / 2) + AOI_WIDTH / 2),
            int(i.lt_pos[1] + y + (rect_size / 2) + AOI_LENGTH / 2),
        )
    i.s_aoi = AOI(
        lt_pos=_lt_pos,
        rt_pos=_rb_pos,
        img=orig_img[_lt_pos[1] : _rb_pos[1], _lt_pos[0] : _rb_pos[0]],
    )

    # 水平方向超出
    if (i.lt_pos[0] + x - (rect_size / 2) - AOI_LENGTH / 2) > i.lt_pos[0]:
        _lt_pos = (
            int(i.lt_pos[0] + x - (rect_size / 2) - AOI_LENGTH / 2),
            int(i.lt_pos[1] + y + (delta / 2) - AOI_WIDTH / 2),
        )
        _rb_pos = (
            int(i.lt_pos[0] + x - (rect_size / 2) + AOI_LENGTH / 2),
            int(i.lt_pos[1] + y + (delta / 2) + AOI_WIDTH / 2),
        )
    else:
        _lt_pos = (
            int(i.lt_pos[0] + x + (rect_size / 2) - AOI_LENGTH / 2),
            int(i.lt_pos[1] + y - (delta / 2) - AOI_WIDTH / 2),
        )
        _rb_pos = (
            int(i.lt_pos[0] + x + (rect_size / 2) + AOI_LENGTH / 2),
            int(i.lt_pos[1] + y - (delta / 2) + AOI_WIDTH / 2),
        )
    i.t_aoi = AOI(
        lt_pos=_lt_pos,
        rt_pos=_rb_pos,
        img=orig_img[_lt_pos[1] : _rb_pos[1], _lt_pos[0] : _rb_pos[0]],
    )

for i in rois:
    if TYPE_CHECKING:
        assert i.t_aoi and i.s_aoi
    cv2.rectangle(img, i.lt_pos, i.rt_pos, (0, 0, 255), 2)
    cv2.rectangle(
        img,
        i.t_aoi.lt_pos,
        i.t_aoi.rt_pos,
        (0, 255, 0),
        2,
    )
    cv2.rectangle(
        img,
        i.s_aoi.lt_pos,
        i.s_aoi.rt_pos,
        (0, 255, 0),
        2,
    )

cv2.namedWindow("result", 0)
cv2.resizeWindow("result", 1600, 1115)  # 设置窗口大小
cv2.imshow('result', img)
# cv2.imwrite('test.png', rois[0].img)
cv2.waitKey()

# if TYPE_CHECKING:
#     assert rois[0].t_aoi and rois[0].s_aoi

# sfr = SFR(verbose=True)
# mtf, status = sfr.calc_sfr(cv2.cvtColor(rois[0].t_aoi.img, cv2.COLOR_BGR2GRAY))

# for i in mtf:
#     print(f'{i[0] / 0.0024}\t:{i[1]}')
# print(status)
